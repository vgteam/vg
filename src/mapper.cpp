#include <unordered_set>
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
    , cluster_min(1)
    , hit_max(0)
    , hit_size_threshold(512)
    , kmer_min(0)
    , kmer_sensitivity_step(5)
    , thread_extension(10)
    , max_thread_gap(30)
    , context_depth(1)
    , max_multimaps(1)
    , max_attempts(0)
    , softclip_threshold(0)
    , max_softclip_iterations(10)
    , prefer_forward(false)
    , greedy_accept(false)
    , accept_identity(0.75)
    , min_identity(0)
    , min_kmer_entropy(0)
    , debug(false)
    , alignment_threads(1)
    , min_mem_length(8)
    , mem_threading(false)
    , max_target_factor(128)
    , max_query_graph_ratio(128)
    , extra_multimaps(1)
    , always_rescue(false)
    , fragment_size(0)
    , fragment_max(1e5)
    , fragment_sigma(4)
    , fragment_length_cache_size(1000)
    , cached_fragment_length_mean(0)
    , cached_fragment_length_stdev(0)
    , since_last_fragment_length_estimate(0)
    , fragment_length_estimate_interval(100)
    , perfect_pair_identity_threshold(0.95)
    , mapping_quality_method(Approx)
    , adjust_alignments_for_base_quality(false)
    , full_length_alignment_bonus(5)
    , max_mapping_quality(64)
    , mem_reseed_length(64)
{
    init_aligner(default_match, default_mismatch, default_gap_open, default_gap_extension);
    init_node_cache();
    init_node_pos_cache();
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
    for (auto& aligner : qual_adj_aligners) {
        delete aligner;
    }
    for (auto& aligner : regular_aligners) {
        delete aligner;
    }
    for (auto& nc : node_cache) {
        delete nc;
    }
    for (auto& np : node_pos_cache) {
        delete np;
    }
}
    
double Mapper::estimate_gc_content() {
    
    uint64_t at = 0, gc = 0;
    
    if (gcsa) {
        at = gcsa::Range::length(gcsa->find(string("A"))) + gcsa::Range::length(gcsa->find(string("T")));
        gc = gcsa::Range::length(gcsa->find(string("G"))) + gcsa::Range::length(gcsa->find(string("C")));
    }
    else if (index) {
        at = index->approx_size_of_kmer_matches("A") + index->approx_size_of_kmer_matches("T");
        gc = index->approx_size_of_kmer_matches("G") + index->approx_size_of_kmer_matches("C");
    }

    if (at == 0 || gc == 0) {
        return default_gc_content;
    }
    
    return ((double) gc) / (at + gc);
}

void Mapper::init_node_cache(void) {
    for (auto& nc : node_cache) {
        delete nc;
    }
    node_cache.clear();
    for (int i = 0; i < alignment_threads; ++i) {
        node_cache.push_back(new LRUCache<id_t, Node>(256));
    }
}

void Mapper::init_node_pos_cache(void) {
    for (auto& nc : node_pos_cache) {
        delete nc;
    }
    node_pos_cache.clear();
    for (int i = 0; i < alignment_threads; ++i) {
        node_pos_cache.push_back(new LRUCache<gcsa::node_type, map<string, vector<size_t> > >(256));
    }
}

void Mapper::clear_aligners(void) {
    for (auto& aligner : qual_adj_aligners) {
        delete aligner;
    }
    qual_adj_aligners.clear();
    for (auto& aligner : regular_aligners) {
        delete aligner;
    }
    regular_aligners.clear();
}

void Mapper::init_aligner(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend) {
    // hacky, find max score so that scaling doesn't change score
    int32_t max_score = match;
    if (mismatch > max_score) max_score = mismatch;
    if (gap_open > max_score) max_score = gap_open;
    if (gap_extend > max_score) max_score = gap_extend;
    
    double gc_content = estimate_gc_content();

    qual_adj_aligners.resize(alignment_threads);
    regular_aligners.resize(alignment_threads);
    for (int i = 0; i < alignment_threads; ++i) {
        qual_adj_aligners[i] = new QualAdjAligner(match, mismatch, gap_open, gap_extend, max_score,
                                         255, gc_content);
        regular_aligners[i] = new Aligner(match, mismatch, gap_open, gap_extend);
        regular_aligners[i]->init_mapping_quality(gc_content); // should be done in constructor
    }
}

void Mapper::set_alignment_scores(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend) {
    if (!qual_adj_aligners.empty()) {
        auto aligner = regular_aligners.front();
        // we've already set the right score
        if (match == aligner->match && mismatch == aligner->mismatch &&
            gap_open == aligner->gap_open && gap_extend == aligner->gap_extension) {
            return;
        }
        // otherwise, destroy them and reset
        clear_aligners();
    }
    // reset the aligners
    init_aligner(match, mismatch, gap_open, gap_extend);
}

// todo add options for aligned global and pinned
Alignment Mapper::align_to_graph(const Alignment& aln,
                                 VG& vg,
                                 size_t max_query_graph_ratio,
                                 bool pinned_alignment,
                                 bool pin_left,
                                 int8_t full_length_bonus,
                                 bool banded_global) {
    // check if we have a cached aligner for this thread
    if (aln.quality().empty()) {
        auto aligner = get_regular_aligner();
        //aligner.align_global_banded(aln, graph.graph, band_padding);
        return vg.align(aln,
                        aligner,
                        max_query_graph_ratio,
                        pinned_alignment,
                        pin_left,
                        full_length_bonus,
                        banded_global);
    } else {
        auto aligner = get_qual_adj_aligner();
        if (adjust_alignments_for_base_quality) {
            return vg.align_qual_adjusted(aln,
                                          aligner,
                                          max_query_graph_ratio,
                                          pinned_alignment,
                                          pin_left,
                                          full_length_bonus,
                                          banded_global);
        } else {
            return vg.align(aln,
                            aligner,
                            max_query_graph_ratio,
                            pinned_alignment,
                            pin_left,
                            full_length_bonus,
                            banded_global);
        }
    }
}

Alignment Mapper::align(const string& seq, int kmer_size, int stride, int max_mem_length, int band_width) {
    Alignment aln;
    aln.set_sequence(seq);
    return align(aln, kmer_size, stride, max_mem_length, band_width);
}

// align read2 near read1's mapping location
void Mapper::align_mate_in_window(const Alignment& read1, Alignment& read2, int pair_window) {
    if (read1.score() == 0) return; // bail out if we haven't aligned the first
    // try to recover in region
    auto& path = read1.path();
    int64_t idf = path.mapping(0).position().node_id();
    int64_t idl = path.mapping(path.mapping_size()-1).position().node_id();
    if(idf > idl) {
        swap(idf, idl);
    }
    // but which way should we expand? this will make things much easier.
    
    // We'll look near the leftmost and rightmost nodes, but we won't try and
    // bridge the whole area of the read, because there may be an ID
    // discontinuity.
    int64_t first = max((int64_t)0, idf - pair_window);
    int64_t last = idl + (int64_t) pair_window;
    
    // Now make sure the ranges don't overlap, because if they do we'll
    // duplicate nodes.
    
    // They can only overlap as idf on top of idl, since we swapped them above.
    // TODO: account at all for orientation? Maybe our left end is in higher
    // numbers but really does need to look left and not right.
    if(idf >= idl) {
        idf--;
    }
    
    VG* graph = new VG;

    if(debug) {
        cerr << "Rescuing in " << first << "-" << idf << " and " << idl << "-" << last << endl;
    }
    
    // TODO: how do we account for orientation when using ID ranges?

    // Now we need to get the neighborhood by ID and expand outward by actual
    // edges. How we do this depends on what indexing structures we have.
    if(xindex) {
        // should have callback here
        xindex->get_id_range(first, idf, graph->graph);
        xindex->get_id_range(idl, last, graph->graph);
        
        // don't get the paths (this isn't yet threadsafe in sdsl-lite)
        xindex->expand_context(graph->graph, context_depth, false);
        graph->rebuild_indexes();
    } else if(index) {
        index->get_range(first, idf, *graph);
        index->get_range(idl, last, *graph);
        index->expand_context(*graph, context_depth);
    } else {
        cerr << "error:[vg::Mapper] cannot align mate with no graph data" << endl;
        exit(1);
    }


    graph->remove_orphan_edges();
    
    if(debug) {
        cerr << "Rescue graph size: " << graph->size() << endl;
    }
    
    read2.clear_path();
    read2.set_score(0);

    read2 = align_to_graph(read2, *graph, max_query_graph_ratio);
    delete graph;
}

map<string, double> Mapper::alignment_mean_path_positions(const Alignment& aln, bool first_hit_only) {
    map<string, double> mean_pos;
    // Alignments are consistent if their median node id positions are within the fragment_size
    
    // We need the sets of nodes visited by each alignment
    set<id_t> ids;
    
    for(size_t i = 0; i < aln.path().mapping_size(); i++) {
        // Collect all the unique nodes visited by the first algnment
        ids.insert(aln.path().mapping(i).position().node_id());
    }
    map<string, map<int, vector<id_t> > > node_positions;
    for(auto id : ids) {
        for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0))) {
            auto& name = ref.first;
            for (auto pos : ref.second) {
                node_positions[name][pos].push_back(id);
            }
        }
        // just get the first one
        if (first_hit_only && node_positions.size()) break;
    }
    // get median mapping positions
    int idscount = 0;
    double idssum = 0;
    for (auto& ref : node_positions) {
        for (auto& p : ref.second) {
            for (auto& n : p.second) {
                auto pos = p.first + get_node_length(n)/2;
                if (ids.count(n)) {
                    idscount++;
                    idssum += pos;
                }
            }
        }
        mean_pos[ref.first] = idssum/idscount;
    }
    return mean_pos;
}

bool Mapper::alignments_consistent(const map<string, double>& pos1,
                                   const map<string, double>& pos2,
                                   int fragment_size_bound) {
    set<string> comm_refs;
    for (auto& p : pos1) {
        auto& name = p.first;
        if (pos2.find(name) != pos2.end()) {
            comm_refs.insert(name);
        }
    }
    // Alignments are consistent if their median node id positions are within the fragment_size
    
    // get median mapping positions
    for (auto& ref : comm_refs) {
        // this is unsafe looking, but we know they have the same keys for these values
        auto mean1 = pos1.find(ref)->second;
        auto mean2 = pos2.find(ref)->second;
        if (abs(mean1 - mean2) < fragment_size_bound) {
            return true;
        }
    }
    return false;
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi(
    const Alignment& read1,
    const Alignment& read2,
    bool& queued_resolve_later,
    int kmer_size,
    int stride,
    int max_mem_length,
    int band_width,
    int pair_window,
    bool only_top_scoring_pair) {

    // use mem threading if requested and we have not need to band (not implemented)
    if (mem_threading && read1.sequence().size() < band_width) {
        return align_paired_multi_simul(read1,
                                        read2,
                                        queued_resolve_later,
                                        max_mem_length,
                                        only_top_scoring_pair);
    } else {
        return align_paired_multi_sep(read1,
                                      read2,
                                      queued_resolve_later,
                                      kmer_size,
                                      stride,
                                      max_mem_length,
                                      band_width,
                                      pair_window,
                                      only_top_scoring_pair);
    }
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi_sep(
    const Alignment& read1,
    const Alignment& read2,
    bool& queued_resolve_later,
    int kmer_size,
    int stride,
    int max_mem_length,
    int band_width,
    int pair_window,
    bool only_top_scoring_pair) {
    // We have some logic around align_mate_in_window to handle orientation
    // Since we now support reversing edges, we have to at least try opposing orientations for the reads.
    auto align_mate = [&](const Alignment& read, Alignment& mate) {
        // Make an alignment to align in the same local orientation as the read
        Alignment aln_same = mate;
        aln_same.clear_path();
        // And one to align in the opposite local orientation
        // Always reverse the opposite direction sequence
        Alignment aln_opposite = reverse_complement_alignment(aln_same, [&](id_t id) {return get_node_length(id);});

        // We can't rescue off an unmapped read
        assert(read.has_path() && read.path().mapping_size() > 0);

        // Do both the alignments
        align_mate_in_window(read, aln_same, pair_window);
        align_mate_in_window(read, aln_opposite, pair_window);

        if(aln_same.score() >= aln_opposite.score()) {
            // TODO: we should prefer opposign local orientations, but we can't
            // really measure them well.
            mate = aln_same;
        } else {
            // Flip the winning reverse alignment back to the original read orientation
            aln_opposite = reverse_complement_alignment(aln_opposite, [&](id_t id) {
                    return get_node_length(id);
                });
            mate = aln_opposite;
        }
    };

    // find the MEMs for the alignments
    vector<MaximalExactMatch> mems1 = find_mems(read1.sequence().begin(),
                                                read1.sequence().end(),
                                                max_mem_length,
                                                mem_reseed_length);
    for (auto& mem : mems1) { get_mem_hits_if_under_max(mem); }
    vector<MaximalExactMatch> mems2 = find_mems(read2.sequence().begin(),
                                                read2.sequence().end(),
                                                max_mem_length,
                                                mem_reseed_length);
    for (auto& mem : mems2) { get_mem_hits_if_under_max(mem); }
    //cerr << "mems before " << mems1.size() << " " << mems2.size() << endl;
    // Do the initial alignments, making sure to get some extras if we're going to check consistency.

    vector<MaximalExactMatch> pairable_mems1, pairable_mems2;
    vector<MaximalExactMatch>* pairable_mems_ptr_1 = nullptr;
    vector<MaximalExactMatch>* pairable_mems_ptr_2 = nullptr;

    // sensitivity ramp
    // first try to get a consistent pair
    // if none is found, re-run the MEM generation with a shorter MEM length
    // if still none is found, align independently with full MEMS; report inconsistent pair
    // optionally use local resolution ...

    // wishlist
    // break out the entire MEM determination logic
    // and merge it with the clustering
    //

    // find the MEMs for the alignments
    if (fragment_size) {
        // use pair resolution filterings on the SMEMs to constrain the candidates
        set<MaximalExactMatch*> pairable_mems = resolve_paired_mems(mems1, mems2);
        for (auto& mem : mems1) if (pairable_mems.count(&mem)) pairable_mems1.push_back(mem);
        for (auto& mem : mems2) if (pairable_mems.count(&mem)) pairable_mems2.push_back(mem);
        pairable_mems_ptr_1 = &pairable_mems1;
        pairable_mems_ptr_2 = &pairable_mems2;
    } else {
        pairable_mems_ptr_1 = &mems1;
        pairable_mems_ptr_2 = &mems2;
    }

    //cerr << pairable_mems1.size() << " and " << pairable_mems2.size() << endl;

    bool report_consistent_pairs = (bool) fragment_size;

    // use MEM alignment on the MEMs matching our constraints
    // We maintain the invariant that these two vectors of alignments are sorted
    // by score, descending, as returned from align_multi_internal.
    double cluster_mq1, cluster_mq2; // XXX not enabled
    vector<Alignment> alignments1 = align_multi_internal(false, read1, kmer_size, stride, max_mem_length,
                                                         band_width, cluster_mq1, extra_multimaps, pairable_mems_ptr_1);
    vector<Alignment> alignments2 = align_multi_internal(false, read2, kmer_size, stride, max_mem_length,
                                                         band_width, cluster_mq2, extra_multimaps, pairable_mems_ptr_2);

    size_t best_score1 = 0;
    size_t best_score2 = 0;
    // A nonzero best score means we have any valid alignments of that read.
    for (auto& aln : alignments1) best_score1 = max(best_score1, (size_t)aln.score());
    for (auto& aln : alignments2) best_score2 = max(best_score2, (size_t)aln.score());

    //bool rescue = !mem_threading && fragment_size != 0; // don't try to rescue if we have a defined fragment size
    bool rescue = fragment_size != 0;

    // Rescue only if the top alignment on one side has no mappings
    if(rescue && best_score1 == 0 && best_score2 != 0) {
        // Must rescue 1 off of 2
#ifdef debug_mapper
        if (debug) cerr << "Rescue read 1 off of read 2" << endl;
#endif
        alignments1.clear();

        // We use this to deduplicate rescue alignments based on their
        // serialized Prtotobuf paths. Relies on protobuf serialization being
        // deterministic.
        set<string> found;

        for(auto base : alignments2) {
            if(base.score() == 0 || !base.has_path() || base.path().mapping_size() == 0) {
                // Can't rescue off this
                continue;
            }
            Alignment mate = read1;
            align_mate(base, mate);

            string serialized;
            mate.path().SerializeToString(&serialized);
            if(!found.count(serialized)) {
                // This is a novel alignment
                alignments1.push_back(mate);
                found.insert(serialized);
            }

            if(!always_rescue) {
                // We only want to rescue off the best one, and they're sorted
                break;
            }
        }
    } else if(rescue && best_score1 != 0 && best_score2 == 0) {
        // Must rescue 2 off of 1
#ifdef debug_mapper
        if (debug) cerr << "Rescue read 2 off of read 1" << endl;
#endif
        alignments2.clear();

        // We use this to deduplicate rescue alignments based on their
        // serialized Prtotobuf paths. Relies on protobuf serialization being
        // deterministic.
        set<string> found;

        for(auto base : alignments1) {
            if(base.score() == 0 || !base.has_path() || base.path().mapping_size() == 0) {
                // Can't rescue off this
                continue;
            }
            Alignment mate = read2;
            align_mate(base, mate);

            string serialized;
            mate.path().SerializeToString(&serialized);
            if(!found.count(serialized)) {
                // This is a novel alignment
                alignments2.push_back(mate);
                found.insert(serialized);
            }

            if(!always_rescue) {
                // We only want to rescue off the best one, and they're sorted
                break;
            }
        }
    } else if(always_rescue) {
        // Try rescuing each off all of the other.
        // We need to be concerned about introducing duplicates.

        // We need temp places to hold the extra alignments we make so as to not
        // rescue off rescues.
        vector<Alignment> extra1;
        vector<Alignment> extra2;

        // We use these to deduplicate alignments based on their serialized
        // Prtotobuf paths. Relies of protobuf serialization being
        // deterministic.
        set<string> found1;
        set<string> found2;

        // Fill in the known alignments
        for(auto existing : alignments1) {
            // Serialize each alignment's Path and put the result in the set
            string serialized;
            existing.path().SerializeToString(&serialized);
            found1.insert(serialized);
        }
        for(auto existing : alignments2) {
            // Serialize each alignment's Path and put the result in the set
            string serialized;
            existing.path().SerializeToString(&serialized);
            found2.insert(serialized);
        }

        for(auto base : alignments1) {
            // Do 2 off of 1
            if(base.score() == 0 || !base.has_path() || base.path().mapping_size() == 0) {
                // Can't rescue off this
                continue;
            }
            Alignment mate = read2;
            align_mate(base, mate);

            string serialized;
            mate.path().SerializeToString(&serialized);
            if(!found2.count(serialized)) {
                // This is a novel alignment
                extra2.push_back(mate);
                found2.insert(serialized);
            }
        }

        for(auto base : alignments2) {
            // Do 1 off of 2
            if(base.score() == 0 || !base.has_path() || base.path().mapping_size() == 0) {
                // Can't rescue off this
                continue;
            }
            Alignment mate = read1;
            align_mate(base, mate);

            string serialized;
            mate.path().SerializeToString(&serialized);
            if(!found1.count(serialized)) {
                // This is a novel alignment
                extra1.push_back(mate);
                found1.insert(serialized);
            }
        }

        // Copy over the new unique alignments
        alignments1.insert(alignments1.end(), extra1.begin(), extra1.end());
        alignments2.insert(alignments2.end(), extra2.begin(), extra2.end());
    }

    // Fix up the sorting by score, descending, in case rescues came out
    // better than normal alignments.
    sort(alignments1.begin(), alignments1.end(), [](const Alignment& a, const Alignment& b) {
            return a.score() > b.score();
        });
    sort(alignments2.begin(), alignments2.end(), [](const Alignment& a, const Alignment& b) {
            return a.score() > b.score();
        });

#ifdef debug_mapper
    if (debug) cerr << alignments1.size() << " alignments for read 1, " << alignments2.size() << " for read 2" << endl;
#endif

    pair<vector<Alignment>, vector<Alignment>> results;

    bool found_consistent = false;


    if (fragment_size) {

        map<Alignment*, map<string, double> > aln_pos;
        for (auto& aln : alignments1) {
            aln_pos[&aln] = alignment_mean_path_positions(aln);
        }
        for (auto& aln : alignments2) {
            aln_pos[&aln] = alignment_mean_path_positions(aln);
        }

        // Now we want to emit consistent pairs, in order of decreasing total score.

        // compare pairs by the sum of their individual scores
        // We need this functor thing to make the priority queue work.
        struct ComparePairedAlignmentScores {
            vector<Alignment>& alns_1;
            vector<Alignment>& alns_2;

        public:
            ComparePairedAlignmentScores(vector<Alignment>& alns_1, vector<Alignment>& alns_2) : alns_1(alns_1), alns_2(alns_2) {}
            bool operator()(const pair<int, int> idxs1, const pair<int, int> idxs2) {
                return (alns_1[idxs1.first].score() + alns_2[idxs1.second].score()
                        < alns_1[idxs2.first].score() + alns_2[idxs2.second].score());
            }
        };

        ComparePairedAlignmentScores compare_paired_alignment_scores = ComparePairedAlignmentScores(alignments1, alignments2);

        // think about the pairs being laid out on a grid over the individual end multimaps, sorted in each dimension by score
        // navigate from top left corner outward to add consistent pairs in decreasing score order
        priority_queue<pair<int, int>, vector<pair<int, int>>, ComparePairedAlignmentScores> pair_queue(compare_paired_alignment_scores);
        // keep track of which indices have been checked to avoid checking them twice when navigating from above and from the left
        std::unordered_set<pair<int, int>> considered_pairs;

        pair<vector<Alignment>, vector<Alignment>> consistent_pairs;
        // ensure that there is always an additional pair to compute a mapping quality against
        int num_pairs = max_multimaps >= 2 ? max_multimaps : 2;

        pair_queue.push(make_pair(0, 0));
        while (!pair_queue.empty() && consistent_pairs.first.size() < num_pairs) {
            // get index of remaining pair with highest combined score
            pair<int, int> aln_pair = pair_queue.top();
            pair_queue.pop();


            if (alignments_consistent(aln_pos[&alignments1[aln_pair.first]], aln_pos[&alignments2[aln_pair.second]], fragment_size)) {
                found_consistent = true;
                consistent_pairs.first.push_back(alignments1[aln_pair.first]);
                consistent_pairs.second.push_back(alignments2[aln_pair.second]);

                if(debug) {
                    cerr << "Found consistent pair " << aln_pair.first << ", " << aln_pair.second
                         << " with scores " << alignments1[aln_pair.first].score()
                         << ", " << alignments2[aln_pair.second].score() << endl;
                }

            }

            // add in the two adjacent indices if we haven't already
            pair<int,int> next_aln_pair_down = make_pair(aln_pair.first + 1, aln_pair.second);
            pair<int,int> next_aln_pair_right = make_pair(aln_pair.first, aln_pair.second + 1);
            if (next_aln_pair_down.first < alignments1.size() && considered_pairs.find(next_aln_pair_down) == considered_pairs.end()) {
                pair_queue.push(next_aln_pair_down);
                considered_pairs.insert(next_aln_pair_down);
            }
            if (next_aln_pair_right.second < alignments2.size() && considered_pairs.find(next_aln_pair_right) == considered_pairs.end()) {
                pair_queue.push(next_aln_pair_right);
                considered_pairs.insert(next_aln_pair_right);
            }
        }

        compute_mapping_qualities(consistent_pairs, cluster_mq1+cluster_mq2);

        // remove the extra pair used to compute mapping quality if necessary
        if (consistent_pairs.first.size() > max_multimaps) {
            consistent_pairs.first.resize(max_multimaps);
            consistent_pairs.second.resize(max_multimaps);
        }

        // mark primary and secondary
        for (int i = 0; i < consistent_pairs.first.size(); i++) {
            consistent_pairs.first[i].mutable_fragment_next()->set_name(read2.name());
            consistent_pairs.first[i].set_is_secondary(i > 0);
            consistent_pairs.second[i].mutable_fragment_prev()->set_name(read1.name());
            consistent_pairs.second[i].set_is_secondary(i > 0);
        }

        // zap everything unless the primary alignments are individually top-scoring
        if (only_top_scoring_pair && consistent_pairs.first.size() &&
            (consistent_pairs.first[0].score() < alignments1[0].score() ||
             consistent_pairs.second[0].score() < alignments2[0].score())) {
            consistent_pairs.first.clear();
            consistent_pairs.second.clear();
        }

        if (!consistent_pairs.first.empty()) {
            results = consistent_pairs;
        } else {
            // no consistent pairs found
            // if we can decrease our MEM size
            // clear, to trigger size reduction

            // otherwise, yolo
        }

    } else {

        results = make_pair(alignments1, alignments2);
        compute_mapping_qualities(results, cluster_mq1 + cluster_mq2);

        // Truncate to max multimaps
        if(results.first.size() > max_multimaps) {
            results.first.resize(max_multimaps);
        }
        if(results.second.size() > max_multimaps) {
            results.second.resize(max_multimaps);
        }

        // mark primary and secondary
        for (int i = 0; i < results.first.size(); i++) {
            results.first[i].mutable_fragment_next()->set_name(read2.name());
            results.first[i].set_is_secondary(i > 0);
        }
        for (int i = 0; i < results.second.size(); i++) {
            results.second[i].mutable_fragment_prev()->set_name(read1.name());
            results.second[i].set_is_secondary(i > 0);
        }

    }

    // change the potential set of MEMs by dropping the maximum MEM size
    // this tends to slightly boost sensitivity at minimal cost
    if (results.first.empty()
        || results.second.empty()
        || !results.first.front().score()
        || !results.second.front().score()) {
        //|| fragment_size && !found_consistent) {
        //cerr << "failed alignment" << endl;
        if (kmer_sensitivity_step) {
            int new_mem_max = max((int) min_mem_length,
                                  (int) (max_mem_length ? max_mem_length : gcsa->order()) - kmer_sensitivity_step);
            if (new_mem_max == min_mem_length) {
                // do noting
            } else if (new_mem_max > min_mem_length) {
                //cerr << "trying with " << new_mem_max << endl;
                return align_paired_multi_sep(read1, read2,
                                              queued_resolve_later,
                                              kmer_size, stride,
                                              new_mem_max,
                                              band_width, pair_window);
            }
        }
    }

    // we tried to align
    // if we don't have a fragment_size yet determined
    // and we didn't get a perfect, unambiguous hit on both reads
    // we'll need to try it again later when we do have a fragment_size
    // so store it in a buffer local to this mapper

    // tag the results with their fragment lengths
    // record the lengths in a deque that we use to keep a running estimate of the fragment length distribution
    // we then set the fragment_size cutoff using the moments of the estimated distribution
    bool imperfect_pair = false;
    for (int i = 0; i < min(results.first.size(), results.second.size()); ++i) {
        auto& aln1 = results.first.at(i);
        auto& aln2 = results.second.at(i);
        auto approx_frag_lengths = approx_pair_fragment_length(aln1, aln2);
        for (auto& j : approx_frag_lengths) {
            Path fragment;
            fragment.set_name(j.first);
            fragment.set_length(j.second);
            *aln1.add_fragment() = fragment;
            *aln2.add_fragment() = fragment;
            // if we have a perfect mapping, and we're under our hard fragment length cutoff
            // push the result into our deque of fragment lengths
            if (results.first.size() == 1
                && results.second.size() == 1
                && results.first.front().identity() > perfect_pair_identity_threshold
                && results.second.front().identity() > perfect_pair_identity_threshold
                && j.second < fragment_max) { // hard cutoff
                //cerr << "aln\tperfect alignments" << endl;
                record_fragment_length(j.second);
            } else if (!fragment_size) {
                imperfect_pair = true;
            }
            //cerr << "aln\t" << aln1.name() << "\t" << aln2.name() << "\t" << j.first << "\t" << j.second << "\t"
            //     << cached_fragment_length_mean << "\t" << cached_fragment_length_stdev << endl;
            //<< fragment_length_mean() << "\t" << fragment_length_stdev() << "\t"
        }
    }

    if (imperfect_pair && fragment_max) {
        imperfect_pairs_to_retry.push_back(make_pair(read1, read2));
        results.first.clear();
        results.second.clear();
        // we signal the fact that this isn't a perfect pair, so we don't write it out externally?
        queued_resolve_later = true;
    }

    // remove results we don't need given our requested number of multimaps
    if (results.first.size() > max_multimaps) {
        results.first.resize(max_multimaps);
    }
    if (results.second.size() > max_multimaps) {
        results.second.resize(max_multimaps);
    }

    if(results.first.empty()) {
        results.first.push_back(read1);
        auto& aln = results.first.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }
    if(results.second.empty()) {
        results.second.push_back(read2);
        auto& aln = results.second.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }

    // Make sure to link up alignments even if they aren't mapped.
    for (auto& aln : results.first) {
        aln.mutable_fragment_next()->set_name(read2.name());
    }

    for (auto& aln : results.second) {
        aln.mutable_fragment_prev()->set_name(read1.name());
    }

    return results;
    
}

pair<vector<Alignment>, vector<Alignment>> Mapper::align_paired_multi_simul(
    const Alignment& read1,
    const Alignment& read2,
    bool& queued_resolve_later,
    int max_mem_length,
    bool only_top_scoring_pair) {

#ifdef debug_mapper
    if (debug) cerr << "aligning read 1 " << pb2json(read1) << endl
                    << "aligning read 2 " << pb2json(read2) << endl;
#endif
    
    double avg_node_len = average_node_length();
    int total_multimaps = max_multimaps + extra_multimaps;
    double cluster_mq = 0;

    pair<vector<Alignment>, vector<Alignment>> results;

    // find the MEMs for the alignments
    vector<MaximalExactMatch> mems1 = find_mems(read1.sequence().begin(),
                                                read1.sequence().end(),
                                                max_mem_length,
                                                mem_reseed_length);
    for (auto& mem : mems1) { get_mem_hits_if_under_max(mem); }
#ifdef debug_mapper
    if (debug) cerr << "mems for read 1 " << mems_to_json(mems1) << endl;
#endif
    vector<MaximalExactMatch> mems2 = find_mems(read2.sequence().begin(),
                                                read2.sequence().end(),
                                                max_mem_length,
                                                mem_reseed_length);
    for (auto& mem : mems2) { get_mem_hits_if_under_max(mem); }
#ifdef debug_mapper
    if (debug) cerr << "mems for read 2 " << mems_to_json(mems2) << endl;
#endif

    auto transition_weight = [&](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {

        // set up positions for distance query
        auto& node1 = m1.nodes.front();
        pos_t m1_pos = make_pos_t(
            gcsa::Node::id(node1),
            gcsa::Node::rc(node1),
            gcsa::Node::offset(node1));
        auto& node2 = m2.nodes.front();
        pos_t m2_pos = make_pos_t(
            gcsa::Node::id(node2),
            gcsa::Node::rc(node2),
            gcsa::Node::offset(node2));

        // approximate distance by node lengths
        double approx_distance = (double) abs(id(m1_pos) - id(m2_pos)) * avg_node_len;

        // are the two mems in the same fragment?
        // we handle the distance metric differently in these cases
        if (m1.fragment < m2.fragment) {
            int unique_coverage = m2.length();
            int max_length = (fragment_size ? fragment_size : fragment_max);
            int dist = approx_distance;
#ifdef debug_mapper
            if (debug) cerr << "distance from " << m1_pos << " to " << m2_pos << " = "<< dist << endl;
#endif
            if (dist >= max_length) {
                return -1.0;
            } else {
                if (is_rev(m1_pos) != is_rev(m2_pos)) {
                    m2_pos = reverse(m2_pos, get_node_length(id(m2_pos)));
                }
                // improve on our approximate metric if we have paths to use
                if (xindex->path_count) {
                    dist = abs(xindex->min_distance_in_paths(
                                   id(m1_pos), is_rev(m1_pos), offset(m1_pos),
                                   id(m2_pos), is_rev(m2_pos), offset(m2_pos)));
                }
                if (fragment_size) {
                    return 1.0/(abs(fragment_size - dist)+1) * unique_coverage;
                } else {
                    return 1.0/(abs(fragment_max - dist)+1) * unique_coverage;
                }
            }
        } else if (m1.fragment > m2.fragment) {
            // don't allow going backwards in the threads
            return -1.0;
        } else { //if (m1.fragment == m2.fragment) {
            int max_length = m1.length() + m2.length();
            // find the difference in m1.end and m2.begin
            // find the positional difference in the graph between m1.end and m2.begin
            int unique_coverage = (m1.end < m2.begin ? m1.length() + m2.length() : m2.end - m1.begin);
#ifdef debug_mapper
            if (debug) cerr << "approx distance " << approx_distance << endl;
#endif
            if (approx_distance > max_length) {
                // too far
                return -1.0;
            } else {
                int distance = graph_distance(m1_pos, m2_pos, max_length);
                if (distance == max_length) {
                    return -1.0;
                }
                double jump = (m2.begin - m1.begin) - distance;
                if (jump < 0) {
                    // disable reversings
                    return -1.0;
                } else {
                    if (is_rev(m1_pos) != is_rev(m2_pos)) {
                        // disable inversions
                        return -1.0;
                    } else {
                        // accepted transition
                        return 1.0 / (double)(abs(jump) + 1) * unique_coverage;
                    }
                }
            }
        }

    };

    // build the paired-read MEM markov model
    MEMMarkovModel markov_model({ read1.sequence().size(), read2.sequence().size() }, { mems1, mems2 }, this, transition_weight, 32);
    vector<vector<MaximalExactMatch> > clusters = markov_model.traceback(total_multimaps, debug);

    // now reconstruct the paired fragments from the threads
    // for each thread we accept either both pairs or one fragment or the other
    //
    //  

    auto show_clusters = [&](void) {
        cerr << "clusters: " << endl;
        for (auto& cluster : clusters) {
            cerr << cluster.size() << " MEMs covering " << cluster_coverage(cluster) << " @ ";
            for (auto& mem : cluster) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << ",";
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

#ifdef debug_mapper
    if (debug) {
        cerr << "### clusters:" << endl;
        show_clusters();
    }
#endif

    vector<pair<Alignment, Alignment> > alns;
    //pair<vector<Alignment>, vector<Alignment> > alns;
    int multimaps = 0;
    for (auto& cluster : clusters) {
        if (multimaps > total_multimaps) { break; }
        // break the cluster into two pieces
        vector<MaximalExactMatch> cluster1, cluster2;
        bool seen1=false, seen2=false;
        for (auto& mem : cluster) {
            if (!seen2 && mem.fragment == 1) {
                cluster1.push_back(mem);
                seen1 = true;
            } else if (mem.fragment == 2) {
                cluster2.push_back(mem);
                seen2 = true;
            } else {
                cerr << "vg map error misordered fragments in cluster" << endl;
                assert(false);
            }
        }
        alns.emplace_back();
        auto& p = alns.back();
        if (cluster1.size()) {
            p.first = patch_alignment(mems_to_alignment(read1, cluster1));
        } else {
            p.first = read1;
            p.first.clear_score();
            p.first.clear_identity();
            p.first.clear_path();
        }
        if (cluster2.size()) {
            p.second = patch_alignment(mems_to_alignment(read2, cluster2));
        } else {
            p.second = read2;
            p.second.clear_score();
            p.second.clear_identity();
            p.second.clear_path();
        }

        ++multimaps;
        
#ifdef debug_mapper
    if (debug) { cerr << "patch identity " << patch.identity() << endl; }
#endif
        /*
        if (patch.identity() > min_identity) {
            alns.emplace_back(patch);
            auto& a = alns.back();
            a.set_name(aln.name());
        }
        */
    }
    // sort the aligned pairs
    std::sort(alns.begin(), alns.end(),
              [&](const pair<Alignment, Alignment>& pair1,
                  const pair<Alignment, Alignment>& pair2) {
                  return pair1.first.score() + pair1.second.score()
                      > pair2.first.score() + pair2.second.score();
              });
    
#ifdef debug_mapper
    if (debug) {
        for (auto& p : alns) {
            auto& aln1 = p.first;
            cerr << "cluster aln 1 ------- " << pb2json(aln1) << endl;
            if (!check_alignment(aln1)) {
                cerr << "alignment failure " << pb2json(aln1) << endl;
                assert(false);
            }
            auto& aln2 = p.second;
            cerr << "cluster aln 2 ------- " << pb2json(aln2) << endl;
            if (!check_alignment(aln2)) {
                cerr << "alignment failure " << pb2json(aln2) << endl;
                assert(false);
            }
        }
    }
#endif

    // calculate cluster mapping quality
    cluster_mq = compute_cluster_mapping_quality(clusters, read1.sequence().size() + read2.sequence().size());
#ifdef debug_mapper
    if (debug) cerr << "cluster mq == " << cluster_mq << endl;
#endif

    // rebuild the thing we'll return
    int read1_max_score = 0;
    int read2_max_score = 0;
    for (auto& p : alns) {
        read1_max_score = max(p.first.score(), read1_max_score);
        read2_max_score = max(p.second.score(), read2_max_score);
        results.first.push_back(p.first);
        results.second.push_back(p.second);
    }
    compute_mapping_qualities(results, cluster_mq);

    // remove the extra pair used to compute mapping quality if necessary
    if (results.first.size() > max_multimaps) {
        results.first.resize(max_multimaps);
        results.second.resize(max_multimaps);
    }

    // mark primary and secondary
    for (int i = 0; i < results.first.size(); i++) {
        results.first[i].mutable_fragment_next()->set_name(read2.name());
        results.first[i].set_is_secondary(i > 0);
        results.second[i].mutable_fragment_prev()->set_name(read1.name());
        results.second[i].set_is_secondary(i > 0);
    }

    // optionally zap everything unless the primary alignments are individually top-scoring
    if (only_top_scoring_pair && results.first.size() &&
        (results.first[0].score() < read1_max_score ||
         results.second[0].score() < read2_max_score)) {
        results.first.clear();
        results.second.clear();
    }
    
    // we tried to align
    // if we don't have a fragment_size yet determined
    // and we didn't get a perfect, unambiguous hit on both reads
    // we'll need to try it again later when we do have a fragment_size
    // so store it in a buffer local to this mapper

    // tag the results with their fragment lengths
    // record the lengths in a deque that we use to keep a running estimate of the fragment length distribution
    // we then set the fragment_size cutoff using the moments of the estimated distribution
    bool imperfect_pair = false;
    for (int i = 0; i < min(results.first.size(), results.second.size()); ++i) {
        auto& aln1 = results.first.at(i);
        auto& aln2 = results.second.at(i);
        auto approx_frag_lengths = approx_pair_fragment_length(aln1, aln2);
        for (auto& j : approx_frag_lengths) {
            Path fragment;
            fragment.set_name(j.first);
            fragment.set_length(j.second);
            *aln1.add_fragment() = fragment;
            *aln2.add_fragment() = fragment;
            // if we have a perfect mapping, and we're under our hard fragment length cutoff
            // push the result into our deque of fragment lengths
            if (results.first.size() == 1
                && results.second.size() == 1
                && results.first.front().identity() > perfect_pair_identity_threshold
                && results.second.front().identity() > perfect_pair_identity_threshold
                && j.second < fragment_max) { // hard cutoff
                //cerr << "aln\tperfect alignments" << endl;
                record_fragment_length(j.second);
            } else if (!fragment_size) {
                imperfect_pair = true;
            }
            //cerr << "aln\t" << aln1.name() << "\t" << aln2.name() << "\t" << j.first << "\t" << j.second << "\t"
            //     << cached_fragment_length_mean << "\t" << cached_fragment_length_stdev << endl;
            //<< fragment_length_mean() << "\t" << fragment_length_stdev() << "\t" 
        }
    }

    if (imperfect_pair && fragment_max) {
        imperfect_pairs_to_retry.push_back(make_pair(read1, read2));
        results.first.clear();
        results.second.clear();
        // we signal the fact that this isn't a perfect pair, so we don't write it out externally?
        queued_resolve_later = true;
    }

    if(results.first.empty()) {
        results.first.push_back(read1);
        auto& aln = results.first.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }
    if(results.second.empty()) {
        results.second.push_back(read2);
        auto& aln = results.second.back();
        aln.clear_path();
        aln.clear_score();
        aln.clear_identity();
    }

    // Make sure to link up alignments even if they aren't mapped.
    for (auto& aln : results.first) {
        aln.set_name(read1.name());
        aln.mutable_fragment_next()->set_name(read2.name());
    }

    for (auto& aln : results.second) {
        aln.set_name(read2.name());
        aln.mutable_fragment_prev()->set_name(read1.name());
    }

    return results;

}

// rank the clusters by the number of unique read bases they cover
int Mapper::cluster_coverage(const vector<MaximalExactMatch>& cluster) {
    set<string::const_iterator> seen;
    for (auto& mem : cluster) {
        string::const_iterator c = mem.begin;
        while (c != mem.end) seen.insert(c++);
    }
    return seen.size();
}

double Mapper::compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters,
                                               int read_length) {
    if (clusters.size() == 0) {
        return 0;
    }
    if (clusters.size() == 1) {
        return { (double)max_mapping_quality };
    }
    //size_t read_length = aln.sequence().size();
    // \prod fraction of the read covered / number of hits
    vector<double> weights;
    for (auto& cluster : clusters) {
        //int unique_cov = cluster_coverage(cluster);
        //cerr << "unique cov " << unique_cov << endl;
        // count up the hits
        //int hits = 0;
        weights.emplace_back();
        double& weight = weights.back();
        for (int i = 0; i < cluster.size(); ++i) {
            // for each mem, count half of its coverage with its neighbors towards this metric
            auto& mem = cluster[i];
            int shared_coverage = 0;
            if (i > 0) {
                auto& prev = cluster[i-1];
                if (prev.fragment == mem.fragment) {
                    shared_coverage += (prev.end <= mem.begin ? 0 : prev.end - mem.begin);
                }
            }
            if (i < cluster.size()-1) {
                auto& next = cluster[i+1];
                if (next.fragment == mem.fragment) {
                    shared_coverage += (mem.end <= next.begin ? 0 : mem.end - next.begin);
                }
            }
            //cerr << "Shared coverage " << shared_coverage << " mem.length() " << mem.length() << endl;
            //int unique_coverage = (m1.end < m2.begin ? m1.length() + m2.length() : m2.end - m1.begin);
            weight +=
                (((double)mem.length() - (double)shared_coverage/2)
                 / read_length)
                / mem.match_count;
        }
        //weight = (weight < 0 ? 0 : weight);
        //cerr << "weight " << weight << endl;
    }
    // return the ratio between best and second best as quality
    std::sort(weights.begin(), weights.end(), std::greater<double>());
    if (weights[0] == 0) return 0;
    return min(max_mapping_quality,
               prob_to_phred(weights[1]/weights[0]));
}

double
Mapper::average_node_length(void) {
    return (double) xindex->seq_length / (double) xindex->node_count;
}

// returns the SMEM clusters which are consistent with the distribution of the MEMs in the read
// provided some tolerance +/-
// uses the exact matches to generate as much of each alignment as possible
// then local dynamic programming to fill in the gaps
vector<Alignment>
Mapper::mems_pos_clusters_to_alignments(const Alignment& aln, vector<MaximalExactMatch>& mems, int additional_multimaps, double& cluster_mq) {

    auto aligner = (aln.quality().empty() ? get_regular_aligner() : get_qual_adj_aligner());
    int total_multimaps = max_multimaps + additional_multimaps;
    
    double avg_node_len = average_node_length();
    // go through the ordered single-hit MEMs
    // build the clustering model
    // find the alignments that are the best-scoring walks through it
    auto transition_weight = [&](const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
        // find the difference in m1.end and m2.begin
        // find the positional difference in the graph between m1.end and m2.begin
        int unique_coverage = (m1.end < m2.begin ? m1.length() + m2.length() : m2.end - m1.begin);
        pos_t m1_pos, m2_pos;
        {
            auto& node = m1.nodes.front();
            id_t id = gcsa::Node::id(node);
            size_t offset = gcsa::Node::offset(node);
            bool is_rev = gcsa::Node::rc(node);
            m1_pos = make_pos_t(id, is_rev, offset);
        }
        {
            auto& node = m2.nodes.front();
            id_t id = gcsa::Node::id(node);
            size_t offset = gcsa::Node::offset(node);
            bool is_rev = gcsa::Node::rc(node);
            m2_pos = make_pos_t(id, is_rev, offset);
        }
        // approximate distance by node lengths
        int max_length = m1.length() + m2.length();
        double approx_distance = (double) abs(id(m1_pos) - id(m2_pos)) * avg_node_len;
#ifdef debug_mapper
        if (debug) cerr << "approx distance " << approx_distance << endl;
#endif
        if (approx_distance > max_length) {
            // too far
            return (double)-1.0;
        } else {
            int distance = graph_distance(m1_pos, m2_pos, max_length);
            if (distance == max_length) {
                return (double)-1.0;
            }
            double jump = (m2.begin - m1.begin) - distance;
            if (jump < 0) {
                // disable reversings
                return (double)-1.0;
            } else {
                if (is_rev(m1_pos) != is_rev(m2_pos)) {
                    // disable inversions
                    return (double)-1.0;
                } else {
                    // accepted transition
                    return (double)1.0 / (double)(abs(jump) + 1) * unique_coverage;
                }
            }
        }
    };

    // build the model
    MEMMarkovModel markov_model({ aln.sequence().size() }, { mems }, this, transition_weight, 10);
    vector<vector<MaximalExactMatch> > clusters = markov_model.traceback(total_multimaps, debug);

    auto show_clusters = [&](void) {
        cerr << "clusters: " << endl;
        for (auto& cluster : clusters) {
            cerr << cluster.size() << " MEMs covering " << cluster_coverage(cluster) << " @ ";
            for (auto& mem : cluster) {
                size_t len = mem.begin - mem.end;
                for (auto& node : mem.nodes) {
                    id_t id = gcsa::Node::id(node);
                    size_t offset = gcsa::Node::offset(node);
                    bool is_rev = gcsa::Node::rc(node);
                    cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << "," << mem.fragment << ",";
                    /*
                    for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0, is_rev))) {
                        auto& name = ref.first;
                        for (auto pos : ref.second) {
                            //cerr << name << (is_rev?"-":"+") << pos + offset;
                            cerr << "|" << id << (is_rev ? "-" : "+") << ":" << offset << ",";
                        }
                    }
                    */
                }
                cerr << mem.sequence() << " ";
            }
            cerr << endl;
        }
    };

#ifdef debug_mapper
    if (debug) {
        cerr << "### clusters:" << endl;
        show_clusters();
    }
#endif

    cluster_mq = compute_cluster_mapping_quality(clusters, aln.sequence().size());

    // for up to our required number of multimaps
    // make the perfect-match alignment for the SMEM cluster
    // then fix it up with DP on the little bits between the alignments
    vector<Alignment> alns;
    int multimaps = 0;
    for (auto& cluster : clusters) {
        if (++multimaps > total_multimaps) { break; }
        Alignment partial_alignment = mems_to_alignment(aln, cluster);
        auto patch = patch_alignment(partial_alignment);
#ifdef debug_mapper
        if (debug) { cerr << "patch identity " << patch.identity() << endl; }
#endif
        if (patch.identity() > min_identity) {
            alns.emplace_back(patch);
            auto& a = alns.back();
            a.set_name(aln.name());
        }
    }
#ifdef debug_mapper
    if (debug) {
        for (auto& aln : alns) {
            cerr << "cluster aln ------- " << pb2json(aln) << endl;
        }
        for (auto& aln : alns) {
            if (!check_alignment(aln)) {
                cerr << "alignment failure " << pb2json(aln) << endl;
                assert(false);
            }
        }
    }
#endif
    return alns;
}

VG Mapper::alignment_subgraph(const Alignment& aln, int context_size) {
    set<id_t> nodes;
    auto& path = aln.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        nodes.insert(path.mapping(i).position().node_id());
    }
    VG graph;
    for (auto& node : nodes) {
        *graph.graph.add_node() = xindex->node(node);
    }
    xindex->expand_context(graph.graph, max(1, context_size)); // get connected edges
    graph.rebuild_indexes();
    return graph;
}

// estimate the fragment length as the difference in mean positions of both alignments
map<string, int> Mapper::approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2) {
    map<string, int> lengths;
    auto pos1 = alignment_mean_path_positions(aln1);
    auto pos2 = alignment_mean_path_positions(aln2);
    for (auto& p : pos1) {
        auto x = pos2.find(p.first);
        if (x != pos2.end()) {
            lengths[p.first] = abs(p.second - x->second);
        }
    }
    return lengths;
}

void Mapper::record_fragment_length(int length) {
    fragment_lengths.push_front(length);
    if (fragment_lengths.size() > fragment_length_cache_size) {
        auto last = fragment_lengths.back();
        fragment_lengths.pop_back();
    }
    if (++since_last_fragment_length_estimate > fragment_length_estimate_interval) {
        cached_fragment_length_mean = fragment_length_mean();
        cached_fragment_length_stdev = fragment_length_stdev();
        // set our fragment size cap to the cached mean + 10x the standard deviation
        fragment_size = cached_fragment_length_mean + fragment_sigma * cached_fragment_length_stdev;
        since_last_fragment_length_estimate = 1;
    }
}

double Mapper::fragment_length_stdev(void) {
    return stdev(fragment_lengths);
}

double Mapper::fragment_length_mean(void) {
    double sum = std::accumulate(fragment_lengths.begin(), fragment_lengths.end(), 0.0);
    return sum / fragment_lengths.size();
}

double Mapper::fragment_length_pdf(double length) {
    return normal_pdf(length, cached_fragment_length_mean, cached_fragment_length_stdev);
}

set<MaximalExactMatch*> Mapper::resolve_paired_mems(vector<MaximalExactMatch>& mems1,
                                                    vector<MaximalExactMatch>& mems2) {
    // find the MEMs that are within estimated_fragment_size of each other

    set<MaximalExactMatch*> pairable;

    // do a wide clustering and then do all pairs within each cluster
    // we will use these to determine the alignment strand
    //map<id_t, StrandCounts> node_strands;
    // records a mapping of id->MEMs, for cluster ranking
    map<id_t, vector<MaximalExactMatch*> > id_to_mems;
    // for clustering
    set<id_t> ids1, ids2;
    vector<id_t> ids;

    // run through the mems
    for (auto& mem : mems1) {
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(&mem);
            ids1.insert(id);
            ids.push_back(id);
        }
    }
    for (auto& mem : mems2) {
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(&mem);
            ids2.insert(id);
            ids.push_back(id);
        }
    }
    // remove duplicates
    //std::sort(ids.begin(), ids.end());
    //ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

    // get each hit's path-relative position
    map<string, map<int, vector<id_t> > > node_positions;
    for (auto& id : ids) {
        for (auto& ref : node_positions_in_paths(gcsa::Node::encode(id, 0))) {
            auto& name = ref.first;
            for (auto pos : ref.second) {
                node_positions[name][pos].push_back(id);
            }
        }
    }

    vector<vector<id_t> > clusters;
    for (auto& g : node_positions) {
        //if (g.second.empty()) continue; // should be impossible
        //cerr << g.first << endl;
        clusters.emplace_back();
        int prev = -1;
        for (auto& x : g.second) {
            auto cluster = &clusters.back();
            //auto& prev = clusters.back().back();
            auto curr = x.first;
            if(debug) {
                cerr << "p/c " << prev << " " << curr << endl;
            }
            if (prev != -1) {
                if (curr - prev <= fragment_size) {
                    // in cluster
#ifdef debug_mapper
                    if (debug) {
                        cerr << "in cluster" << endl;
                    }
#endif
                } else {
                    // It's a new cluster
                    clusters.emplace_back();
                    cluster = &clusters.back();
                }
            }
            //cerr << " " << x.first << endl;
            for (auto& y : x.second) {
                //cerr << "  " << y << endl;
                cluster->push_back(y);
            }
            prev = curr;
        }
    }

    for (auto& cluster : clusters) {
        // for each pair of ids in the cluster
        // which are not from the same read
        // estimate the distance between them
        // we're roughly in the expected range
        bool has_first = false;
        bool has_second = false;
        for (auto& id : cluster) {
            has_first |= ids1.count(id);
            has_second |= ids2.count(id);
        }
        if (!has_first || !has_second) continue;
        for (auto& id : cluster) {
            for (auto& memptr : id_to_mems[id]) {
                pairable.insert(memptr);
            }
        }
    }

    return pairable;
}

// We need a function to get the lengths of nodes, in case we need to
// reverse an Alignment, including all its Mappings and Positions. TODO:
// make this cache the node lengths for the nodes used in the actual
// alignments somehow?
int64_t Mapper::get_node_length(int64_t node_id) {
    if(xindex) {
        // Grab the node sequence only from the XG index and get its size.
        return xindex->node_length(node_id);
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

bool Mapper::check_alignment(const Alignment& aln) {
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
        if (aln.quality().size() && aln.quality().size() != aln.sequence().size()) {
            cerr << "alignment quality is not the same length as its sequence" << endl
                 << pb2json(aln) << endl;
            return false;
        }
        if (seq != aln.sequence()) {
            cerr << "alignment does not match graph " << endl
                 << pb2json(aln) << endl
                 << "expect:\t" << aln.sequence() << endl
                 << "got:\t" << seq << endl;
            // save alignment
            write_alignment_to_file(aln, "fail-" + hash_alignment(aln) + ".gam");
            // save graph, bigger fragment
            xindex->expand_context(sub, 5, true);
            VG gn; gn.extend(sub);
            gn.serialize_to_file("fail-" + gn.hash() + ".vg");
            return false;
        }
    }
    return true;
}

Alignment Mapper::align_banded(const Alignment& read, int kmer_size, int stride, int max_mem_length, int band_width) {
    // split the alignment up into overlapping chunks of band_width size
    list<Alignment> alignments;
    // force used bandwidth to be divisible by 4
    // round up so we have > band_width
#ifdef debug_mapper
    if (debug) {
        cerr << "trying band width " << band_width << endl;
    }
#endif
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
#ifdef debug_mapper
    if (debug) {
        cerr << "Segment size be " << segment_size << "/" << read.sequence().size() << endl;
    }
#endif
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

    auto do_band = [&](int i) {
        if (max_multimaps > 1) {
            vector<Alignment>& malns = multi_alns[i];
            double cluster_mq = 0;
            malns = align_multi_internal(false, bands[i], kmer_size, stride, max_mem_length, band_width, cluster_mq, extra_multimaps, nullptr);
            // always include an unaligned mapping
            malns.push_back(bands[i]);
            for (vector<Alignment>::iterator a = malns.begin(); a != malns.end(); ++a) {
                Alignment& aln = *a;
                bool above_threshold = aln.identity() >= min_identity;
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
            aln = align(bands[i], kmer_size, stride, max_mem_length, band_width);
            bool above_threshold = aln.identity() >= min_identity;
            if (!above_threshold) {
                aln = bands[i]; // unmapped
            }
            
#ifdef debug_mapper
            if (debug) {
                cerr << "Unstripped alignment: " << pb2json(aln) << endl;
            }
#endif
            
            // strip overlaps
            //cerr << "checking before strip" << endl;
            //check_alignment(aln);
            aln = strip_from_start(aln, to_strip[i].first);
            aln = strip_from_end(aln, to_strip[i].second);
            //cerr << "checking after strip" << endl;
            //check_alignment(aln);
            //cerr << "OK" << endl;
        }
    };
    
    if (alignment_threads > 1) {
#pragma omp parallel for schedule(dynamic,1)
        for (int i = 0; i < bands.size(); ++i) {
            do_band(i);
        }
    } else {
        for (int i = 0; i < bands.size(); ++i) {
            do_band(i);
        }
    }

    // resolve the highest-scoring traversal of the multi-mappings
    if (max_multimaps > 1) {
        alns = resolve_banded_multi(multi_alns);
        multi_alns.clear(); // clean up
    }

    // check that the alignments are valid
#ifdef debug_mapper
    if (debug) {
        for (auto& aln : alns) {
            check_alignment(aln);
        }
    }
#endif

    // merge the resulting alignments
    Alignment merged = merge_alignments(alns);

    merged.set_score(score_alignment(merged));
    merged.set_identity(identity(merged.path()));
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
#ifdef debug_mapper
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
#endif
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

QualAdjAligner* Mapper::get_qual_adj_aligner(void) {
    int tid = qual_adj_aligners.size() > 1 ? omp_get_thread_num() : 0;
    return qual_adj_aligners[tid];
}

Aligner* Mapper::get_regular_aligner(void) {
    int tid = regular_aligners.size() > 1 ? omp_get_thread_num() : 0;
    return regular_aligners[tid];
}

LRUCache<id_t, Node>& Mapper::get_node_cache(void) {
    int tid = node_cache.size() > 1 ? omp_get_thread_num() : 0;
    return *node_cache[tid];
}

LRUCache<gcsa::node_type, map<string, vector<size_t> > >& Mapper::get_node_pos_cache(void) {
    int tid = node_pos_cache.size() > 1 ? omp_get_thread_num() : 0;
    return *node_pos_cache[tid];
}

void Mapper::compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq) {
    if (alns.empty()) return;
    auto aligner = (alns.front().quality().empty() ? get_regular_aligner() : get_qual_adj_aligner());
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_mapping_quality(alns, max_mapping_quality, cluster_mq, true);
            break;
        case Exact:
            aligner->compute_mapping_quality(alns, max_mapping_quality, cluster_mq, false);
            break;
        default: // None
            break;
    }
}
    
void Mapper::compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq) {
    if (pair_alns.first.empty() || pair_alns.second.empty()) return;
    auto aligner = (pair_alns.first.front().quality().empty() ? get_regular_aligner() : get_qual_adj_aligner());
    switch (mapping_quality_method) {
        case Approx:
            aligner->compute_paired_mapping_quality(pair_alns, max_mapping_quality, cluster_mq, true);
            break;
        case Exact:
            aligner->compute_paired_mapping_quality(pair_alns, max_mapping_quality, cluster_mq, false);
            break;
        default: // None
            break;
    }
}

vector<Alignment> Mapper::score_sort_and_deduplicate_alignments(vector<Alignment>& all_alns, const Alignment& original_alignment) {
    if (all_alns.size() == 0) {
        all_alns.emplace_back();
        Alignment& aln = all_alns.back();
        aln = original_alignment;
        aln.clear_path();
        aln.set_score(0);
        return all_alns;
    }
    
    map<int, set<Alignment*> > alignment_by_score;
    for (auto& ta : all_alns) {
        Alignment* aln = &ta;
        alignment_by_score[aln->score()].insert(aln);
    }
    // TODO: Filter down subject to a minimum score per base or something?
    // Collect all the unique alignments (to compute mapping quality) and order by score
    vector<Alignment> sorted_unique_alignments;
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
                sorted_unique_alignments.push_back(*pointer);
                
                // Save it so we can avoid putting it in the vector again
                serializedAlignmentsUsed.insert(serialized);
            }
        }
    }
    return sorted_unique_alignments;
}
    
// filters down to requested number of alignments and marks
void Mapper::filter_and_process_multimaps(vector<Alignment>& sorted_unique_alignments, int additional_multimaps) {
    int total_multimaps = max_multimaps + additional_multimaps;
    if (sorted_unique_alignments.size() > total_multimaps){
        sorted_unique_alignments.resize(total_multimaps);
    }
    
    // TODO log best alignment score?
    for(size_t i = 0; i < sorted_unique_alignments.size(); i++) {
        // Mark all but the first, best alignment as secondary
        sorted_unique_alignments[i].set_is_secondary(i > 0);
    }
}
    
vector<Alignment> Mapper::align_multi(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width) {
    double cluster_mq = 0;
    return align_multi_internal(true, aln, kmer_size, stride, max_mem_length, band_width, cluster_mq, extra_multimaps, nullptr);
}
    
vector<Alignment> Mapper::align_multi_internal(bool compute_unpaired_quality,
                                               const Alignment& aln,
                                               int kmer_size, int stride,
                                               int max_mem_length,
                                               int band_width,
                                               double& cluster_mq,
                                               int additional_multimaps,
                                               vector<MaximalExactMatch>* restricted_mems) {
    
    if(debug) {
        cerr << "align_multi_internal("
            << compute_unpaired_quality << ", " 
            << aln.sequence() << ", " 
            << kmer_size << ", " 
            << stride << ", " 
            << band_width << ", " 
            << additional_multimaps << ", " 
            << restricted_mems << ")" 
            << endl;
    }
    
    // trigger a banded alignment if we need to
    // note that this will in turn call align_multi_internal on fragments of the read
    if (aln.sequence().size() > band_width) {
        // TODO: banded alignment currently doesn't support mapping qualities because it only produces one alignment
#ifdef debug_mapper
        if (debug) cerr << "switching to banded alignment" << endl;
#endif
        return vector<Alignment>{align_banded(aln, kmer_size, stride, max_mem_length, band_width)};
    }
    
    // try to get at least 2 multimaps so that we can calculate mapping quality
    int additional_multimaps_for_quality;
    if (additional_multimaps == 0 && max_multimaps == 1 && mapping_quality_method != None) {
        additional_multimaps_for_quality = 1;
    }
    else {
        additional_multimaps_for_quality = additional_multimaps;
    }

    vector<Alignment> alignments;
    if (kmer_size || xindex == nullptr) {
        // if we've defined a kmer size, use the legacy style mapper
        alignments = align_multi_kmers(aln, kmer_size, stride, band_width);
    }
    else {
        // otherwise use the mem mapper, which is a banded multi mapper by default
        
        // use pre-restricted mems for paired mapping or find mems here
        if (restricted_mems != nullptr) {
            // mem hits will already have been queried
            alignments = align_mem_multi(aln, *restricted_mems, cluster_mq, additional_multimaps_for_quality);
        }
        else {
            vector<MaximalExactMatch> mems = find_mems(aln.sequence().begin(), aln.sequence().end(), max_mem_length, mem_reseed_length);
            // query mem hits
#ifdef debug_mapper
            if (debug) cerr << "mems before filtering " << mems_to_json(mems) << endl;
#endif
            for (auto& mem : mems) { get_mem_hits_if_under_max(mem); }
#ifdef debug_mapper
            if (debug) cerr << "mems after filtering " << mems_to_json(mems) << endl;
#endif
            alignments = align_mem_multi(aln, mems, cluster_mq, additional_multimaps_for_quality);
        }
    }
    
    alignments = score_sort_and_deduplicate_alignments(alignments, aln);
    
    // compute mapping quality before removing extra alignments
    if (compute_unpaired_quality) {
        compute_mapping_qualities(alignments, cluster_mq);
        filter_and_process_multimaps(alignments, 0);
    } else {
        filter_and_process_multimaps(alignments, additional_multimaps);
    }
    
    return alignments;
}

vector<Alignment> Mapper::align_multi_kmers(const Alignment& aln, int kmer_size, int stride, int band_width) {

    std::chrono::time_point<std::chrono::system_clock> start_both, end_both;
#ifdef debug_mapper
    if (debug) start_both = std::chrono::system_clock::now();
#endif
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

#ifdef debug_mapper
    if (debug) cerr << "aligning " << aln.sequence() << endl;
#endif

    // This will hold the best forward alignment (or an alignment with no path and 0 score if no alignment is found).
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
#ifdef debug_mapper
        if (debug) cerr << "realigning with " << kmer_size << " " << stride << endl;
#endif
    };

    int attempt = 0;
    int kmer_count_f = 0;
    int kmer_count_r = 0;

    while (!(best_f.identity() > min_identity
             || best_r.identity() > min_identity)
           && attempt < max_attempts) {

        {
            std::chrono::time_point<std::chrono::system_clock> start, end;
#ifdef debug_mapper
            if (debug) start = std::chrono::system_clock::now();
#endif
            // Go get all the forward alignments, putting the best one in best_f.
            alignments_f = align_threaded(best_f, kmer_count_f, kmer_size, stride, attempt);
#ifdef debug_mapper
            if (debug) {
                end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                cerr << elapsed_seconds.count() << "\t" << "+" << "\t" << best_f.sequence() << endl;
            }
#endif
        }

        if (!(prefer_forward && best_f.identity() >= accept_identity))
        {
            // If we need to look on the reverse strand, do that too.
            std::chrono::time_point<std::chrono::system_clock> start, end;
#ifdef debug_mapper
            if (debug) start = std::chrono::system_clock::now();
#endif
            auto alns =  align_threaded(best_r, kmer_count_r, kmer_size, stride, attempt);
            alignments_r = reverse_complement_alignments(alns,
                                                         (function<int64_t(int64_t)>) ([&](int64_t id) { return get_node_length(id); }));
#ifdef debug_mapper
            if (debug) {
                end = std::chrono::system_clock::now();
                std::chrono::duration<double> elapsed_seconds = end-start;
                cerr << elapsed_seconds.count() << "\t" << "-" << "\t" << best_r.sequence() << endl;
            }
#endif
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

#ifdef debug_mapper
    if (debug) {
        end_both = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end_both-start_both;
        cerr << elapsed_seconds.count() << "\t" << "b" << "\t" << sequence << endl;
    }
#endif
    
    // merge return all alignments found, don't bother sorting here
    vector<Alignment> merged;
    merged.reserve(alignments_f.size() + alignments_r.size());
    for (int i = 0; i < alignments_f.size(); i++) {
        merged.push_back(alignments_f[i]);
    }
    for (int i = 0; i < alignments_r.size(); i++) {
        merged.push_back(alignments_r[i]);
    }
    
    // Return the merged list of good alignments. Does not bother updating the input alignment.
    return merged;
}



Alignment Mapper::align(const Alignment& aln, int kmer_size, int stride, int max_mem_length, int band_width) {
    // TODO computing mapping quality could be inefficient depending on the method chosen
    
    // Do the multi-mapping
    vector<Alignment> best = align_multi(aln, kmer_size, stride, max_mem_length, band_width);

    if(best.size() == 0) {
        // Spit back an alignment that says we failed, but make sure it has the right sequence in it.
        Alignment failed = aln;
        failed.clear_path();
        failed.set_score(0);
        return failed;
    }

    // Otherwise, just report the best alignment, since we know one exists
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
Mapper::find_mems(string::const_iterator seq_begin,
                  string::const_iterator seq_end,
                  int max_mem_length,
                  int reseed_length) {
    
    if (!gcsa) {
        cerr << "error:[vg::Mapper] a GCSA2 index is required to query MEMs" << endl;
        exit(1);
    }

    string::const_iterator cursor = seq_end;
    vector<MaximalExactMatch> mems;

    // an empty sequence matches the entire bwt
    if (seq_begin == seq_end) {
        mems.emplace_back(
            MaximalExactMatch(seq_begin, seq_end,
                              gcsa::range_type(0, gcsa->size() - 1)));
        return mems;
    }
    
    // find SMEMs using GCSA+LCP array
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
    auto full_range = gcsa::range_type(0, gcsa->size() - 1);
    MaximalExactMatch match(cursor, cursor, full_range);
    gcsa::range_type last_range = match.range;
    --cursor; // start off looking at the last character in the query
    while (cursor >= seq_begin) {
        // hold onto our previous range
        last_range = match.range;
        // execute one step of LF mapping
        match.range = gcsa->LF(match.range, gcsa->alpha.char2comp[*cursor]);
        if (gcsa::Range::empty(match.range)
            || max_mem_length && match.end-cursor > max_mem_length
            || match.end-cursor > gcsa->order()) {
            // break on N; which for DNA we assume is non-informative
            // this *will* match many places in assemblies; this isn't helpful
            if (*cursor == 'N' || last_range == full_range) {
                // we mismatched in a single character
                // there is no MEM here
                match.begin = cursor+1;
                match.range = last_range;
                mems.push_back(match);
                match.end = cursor;
                match.range = full_range;
                --cursor;
            } else {
                // we've exhausted our BWT range, so the last match range was maximal
                // or: we have exceeded the order of the graph (FPs if we go further)
                //     we have run over our parameter-defined MEM limit
                // record the last MEM
                match.begin = cursor+1;
                match.range = last_range;
                mems.push_back(match);
                // set up the next MEM using the parent node range
                // length of last MEM, which we use to update our end pointer for the next MEM
                size_t last_mem_length = match.end - match.begin;
                // get the parent suffix tree node corresponding to the parent of the last MEM's STNode
                gcsa::STNode parent = lcp->parent(last_range);
                // change the end for the next mem to reflect our step size
                size_t step_size = last_mem_length - parent.lcp();
                match.end = mems.back().end-step_size;
                // and set up the next MEM using the parent node range
                match.range = parent.range();
            }
        } else {
            // we are matching
            match.begin = cursor;
            // just step to the next position
            --cursor;
        }
    }
    // if we have a non-empty MEM at the end, record it
    if (match.end - match.begin > 0) mems.push_back(match);

    // find the SMEMs from the mostly-SMEM and some MEM list we've built
    // FIXME: un-hack this (it shouldn't be needed!)
    // the algorithm sometimes generates MEMs contained in SMEMs
    // with the pattern that they have the same beginning position
    map<string::const_iterator, string::const_iterator> smems_begin;
    for (auto& mem : mems) {
        auto x = smems_begin.find(mem.begin);
        if (x == smems_begin.end()) {
            smems_begin[mem.begin] = mem.end;
        } else {
            if (x->second < mem.end) {
                x->second = mem.end;
            }
        }
    }
    // remove zero-length entries and MEMs that aren't SMEMs
    // the zero-length ones are associated with single-base MEMs that tend to
    // match the entire index (typically Ns)
    // minor TODO: fix the above algorithm so they aren't introduced at all
    mems.erase(std::remove_if(mems.begin(), mems.end(),
                              [&smems_begin](const MaximalExactMatch& m) {
                                  return m.end-m.begin == 0
                                      || smems_begin[m.begin] != m.end;
                              }),
               mems.end());
    // return the matches in natural order
    std::reverse(mems.begin(), mems.end());
    
    // reseed the long smems with shorter mems
    if (reseed_length) {
        // find if there are any mems that should be reseeded
        // iterate through MEMs
        vector<MaximalExactMatch> reseeded;
        for (auto& mem : mems) {
            if (mem.length() >= reseed_length) {
                // get the count
                mem.fill_match_count(gcsa);
                // only reseed when there is a single hit
                if (mem.match_count > 1) {
                    continue;
                }
                vector<MaximalExactMatch> remems = find_mems(mem.begin, mem.end, mem.length()/2, 0);
                for (auto& rmem : remems) {
                    rmem.fill_match_count(gcsa);
                    // keep if we have more than the match count of the parent
                    if (rmem.match_count > mem.match_count) {
                        reseeded.push_back(rmem);
                    }
                }
            }
        }
        // extend
        mems.reserve(mems.size() + distance(reseeded.begin(), reseeded.end()));
        mems.insert(mems.end(), reseeded.begin(), reseeded.end());
        // re-sort the MEMs by their start position
        std::sort(mems.begin(), mems.end(), [](const MaximalExactMatch& m1, const MaximalExactMatch& m2) { return m1.begin < m2.begin; });
    }
    
    // verify the matches (super costly at scale)
#ifdef debug_mapper
    if (debug) { check_mems(mems); }
#endif
    return mems;
}

void Mapper::check_mems(const vector<MaximalExactMatch>& mems) {
    for (auto mem : mems) {
        get_mem_hits_if_under_max(mem);
        if (sequence_positions(mem.sequence()) != gcsa_nodes_to_positions(mem.nodes)) {
            cerr << "mem failed! " << mem.sequence()
                 << " expected " << sequence_positions(mem.sequence()).size() << " hits "
                 << "but found " << gcsa_nodes_to_positions(mem.nodes).size() << endl;
        }
    }
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

char Mapper::pos_char(pos_t pos) {
    return xg_cached_pos_char(pos, xindex, get_node_cache());
}

map<pos_t, char> Mapper::next_pos_chars(pos_t pos) {
    return xg_cached_next_pos_chars(pos, xindex, get_node_cache());
}

int Mapper::graph_distance(pos_t pos1, pos_t pos2, int maximum) {
    return xg_cached_distance(pos1, pos2, xindex, get_node_cache(), maximum);
}

set<pos_t> Mapper::positions_bp_from(pos_t pos, int distance, bool rev) {
    return xg_cached_positions_bp_from(pos, distance, rev, xindex, get_node_cache());
}

// use LRU caching to get the most-recent node positions
map<string, vector<size_t> > Mapper::node_positions_in_paths(gcsa::node_type node) {
    auto& pos_cache = get_node_pos_cache();
    auto cached = pos_cache.retrieve(node);
    if(!cached.second) {
        // todo use approximate estimate
        cached.first = xindex->position_in_paths(gcsa::Node::id(node), gcsa::Node::rc(node), gcsa::Node::offset(node));
        pos_cache.put(node, cached.first);
    }
    return cached.first;
}

Alignment Mapper::walk_match(const string& seq, pos_t pos) {
    //cerr << "in walk match with " << seq << " " << seq.size() << " " << pos << endl;
    Alignment aln;
    aln.set_sequence(seq);
    auto alns = walk_match(aln, seq, pos);
    if (!alns.size()) {
        //cerr << "no alignments returned from walk match with " << seq << " " << seq.size() << " " << pos << endl;
        //assert(false);
        return aln;
    }
    aln = alns.front(); // take the first one we found
    //assert(alignment_to_length(aln) == alignment_from_length(aln));
    if (alignment_to_length(aln) != alignment_from_length(aln)
        || alignment_to_length(aln) != seq.size()) {
        //cerr << alignment_to_length(aln) << " is not " << seq.size() << endl;
        //cerr << pb2json(aln) << endl;
        //assert(false);
        aln.clear_path();
    }
#ifdef debug_mapper
    if (debug) {
        cerr << "walk_match result " << pb2json(aln) << endl;
        if (!check_alignment(aln)) {
            cerr << "aln is invalid!" << endl;
            exit(1);
        }
    }
#endif
    return aln;
}

vector<Alignment> Mapper::walk_match(const Alignment& base, const string& seq, pos_t pos) {
    //cerr << "in walk_match " << seq << " from " << pos << " with base " << pb2json(base) << endl;
    // go to the position in the xg index
    // and step in the direction given
    // until we exhaust our sequence
    // or hit another node
    vector<Alignment> alns;
    Alignment aln = base;
    Path& path = *aln.mutable_path();
    Mapping* mapping = path.add_mapping();
    *mapping->mutable_position() = make_position(pos);
    // get the first node we match
    int total = 0;
    size_t match_len = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        char c = seq[i];
        //cerr << string(base.path().mapping_size(), ' ') << pos << " @ " << i << " on " << c << endl;
        auto nexts = next_pos_chars(pos);
        // we can have a match on the current node
        if (nexts.size() == 1 && id(nexts.begin()->first) == id(pos)) {
            pos_t npos = nexts.begin()->first;
            // check that the next position would match
            if (i+1 < seq.size()) {
                // we can't step, so we break
                //cerr << "Checking if " << pos_char(npos) << " != " << seq[i+1] << endl;
                if (pos_char(npos) != seq[i+1]) {
                    //cerr << "broke" << endl;
                    //cerr << "returning ..." << alns.size() << endl;
                    return alns;
                }
            }
            // otherwise we step our counters
            ++match_len;
            ++get_offset(pos);
        } else { // or we go into the next node
            // we must be going into another node
            // emit the mapping for this node
            //cerr << "we are going into a new node" << endl;
            // finish the last node
            {
                // we must have matched / we already checked
                ++match_len;
                Edit* edit = mapping->add_edit();
                edit->set_from_length(match_len);
                edit->set_to_length(match_len);
                // reset our counter
                match_len = 0;
            }
            // find the next node that matches our MEM
            bool got_match = false;
            if (i+1 < seq.size()) {
                //cerr << "nexts @ " << i << " " << nexts.size() << endl;
                for (auto& p : nexts) {
                    //cerr << " next : " << p.first << " " << p.second << " (looking for " << seq[i+1] << ")" << endl;
                    if (p.second == seq[i+1]) {
                        if (!got_match) {
                            pos = p.first;
                            got_match = true;
                        } else {
                            auto v = walk_match(aln, seq.substr(i+1), p.first);
                            if (v.size()) {
                                alns.reserve(alns.size() + distance(v.begin(), v.end()));
                                alns.insert(alns.end(), v.begin(), v.end());
                            }
                        }
                    }
                }
                if (!got_match) {
                    // this matching ends here
                    // and we haven't finished matching
                    // thus this path doesn't contain the match
                    //cerr << "got no match" << endl;
                    return alns;
                }

                // set up a new mapping
                mapping = path.add_mapping();
                *mapping->mutable_position() = make_position(pos);
            } else {
                //cerr << "done!" << endl;
            }
        }
    }
    if (match_len) {
        Edit* edit = mapping->add_edit();
        edit->set_from_length(match_len);
        edit->set_to_length(match_len);
    }
    alns.push_back(aln);
    //cerr << "returning " << alns.size() << endl;
    return alns;
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

Alignment Mapper::patch_alignment(const Alignment& aln) {
#ifdef debug_mapper
    if (debug) {
        cerr << "patching " << pb2json(aln) << endl;
        if (!check_alignment(aln)) {
            cerr << "aln is invalid!" << endl;
            exit(1);
        }
    }
#endif
    Alignment patched;
    int score = 0;
    // walk along the alignment and find the portions that are unaligned
    int read_pos = 0;
    auto& path = aln.path();
    auto aligner = aln.quality().empty() ? get_regular_aligner() : get_qual_adj_aligner();
    auto qual_adj_aligner = (QualAdjAligner*) aligner;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
#ifdef debug_mapper
        if (debug) cerr << "looking at mapping " << pb2json(mapping) << endl;
#endif
        pos_t ref_pos = make_pos_t(mapping.position());
        Mapping* new_mapping = patched.mutable_path()->add_mapping();
        *new_mapping->mutable_position() = mapping.position();
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
#ifdef debug_mapper
            if (debug) cerr << "looking at edit " << pb2json(edit) << endl;
#endif
            if (edit_is_match(edit)) {
                // matches behave as expected
                //cerr << "edit is match" << endl;
                if (!aln.quality().empty()) {
                    //cerr << read_pos << " " << edit.to_length() << endl;
                    score += qual_adj_aligner->score_exact_match(
                        aln.sequence().substr(read_pos, edit.to_length()),
                        aln.quality().substr(read_pos, edit.to_length()));
                } else {
                    score += edit.from_length()*aligner->match;
                }
                *new_mapping->add_edit() = edit;
            } else if (edit_is_deletion(edit)) {
                // we can't do anything for deletions-- anyway they shouldn't get here if we call this
                // in the SMEM threading alignment
                score -= aligner->gap_open + edit.from_length()*aligner->gap_extension;
                *new_mapping->add_edit() = edit;
            } else if (edit_is_insertion(edit)) {
                //cerr << "looking at " << edit.sequence() << endl;
                // bits to patch in are recorded like insertions
                // pick up the graph from the start to end where we have an unaligned bit
                // but bail out if we get a lot of graph
                bool go_forward = !is_rev(ref_pos);
                bool go_backward = is_rev(ref_pos);
                id_t id1 = id(ref_pos);
                id_t id2 = 0;
                pos_t after_pos = ref_pos;
                bool soft_clip_to_left = false;
                bool soft_clip_to_right = false;
                // this is a soft clip
                if (i == 0 && j == 0) {
                    //cerr << "first soft clip" << endl;
                    // todo we should flip the orientation of the soft clip flag around if we are reversed
                    // ...
                    //soft_clip_on_start = true;
                    if (is_rev(ref_pos)) {
                        soft_clip_to_right = true;
                        go_forward = true;
                        go_backward = false;
                    } else {
                        soft_clip_to_left = true;
                        go_forward = false;
                        go_backward = true;
                    }
                } else if (j+1 < mapping.edit_size()) {
                    id2 = id1;
                    //cerr << "more edits to go on this node: " << id2 << " " << after_pos << endl;
                } else if (i+1 < path.mapping().size()) {
                    // get the next position in the partial alignment we're patching
                    id2 = path.mapping(i+1).position().node_id();
                    after_pos = make_pos_t(path.mapping(i+1).position());
                    //cerr << "we're up so look for id: " << id2 << " " << after_pos << endl;
                } else {
                    //cerr << "last soft clip" << endl;
                    if (is_rev(ref_pos)) {
                        soft_clip_to_left = true;
                        go_forward = false;
                        go_backward = true;
                    } else {
                        soft_clip_to_right = true;
                        go_forward = true;
                        go_backward = false;
                    }
                }
                //cerr << "working from " << ref_pos << endl;
                // only go backward if we are at the first edit (e.g. soft clip)
                // otherwise we go forward

                // find the cut positions (on the forward strand)
                // if they are on the same node, use the multi-cut interface
                // if they are on separate, cut each node
                // todo... update the multi-cut interface to produce a translation
                // or alternatively, write the translation here for prototyping
                // generate the translation for the rest of the graph
                // delete any head or tail bits that we shouldn't be able to align to
                // instantiate a translator object (make sure it can handle the deleted bits?)
                // translate the alignment

                pos_t first_cut = ref_pos;
                pos_t second_cut = ref_pos;
                bool insertion_between_mems = false;
                if (j+1 < mapping.edit_size()) {
                    //cerr << "not the last edit" << endl;
                    if (edit.from_length()) {
                        get_offset(second_cut) += edit.from_length();
                    } else if (i != 0 && j != 0) {
                        insertion_between_mems = true;
                    }
                } else if (i+1 < path.mapping_size()) {
                    //cerr << "not the last mapping" << endl;
                    // we have to look at the next mapping
                    second_cut = make_pos_t(path.mapping(i+1).position());
                } else {
                    //cerr << "end of alignment" << endl;
                    // nothing to do
                }
                pos_t next_pos = second_cut;

                //cerr << "first_cut before " << first_cut << endl;
                //cerr << "second_cut before " << second_cut << endl;

                // if we get a target graph
                int min_distance = edit.to_length() * 3;

                //cerr << "going at least " << min_distance << endl;
                VG graph;
                if (!insertion_between_mems) {
                    xindex->get_id_range(id1, id1, graph.graph);
                    xindex->expand_context(graph.graph,
                                           min_distance,
                                           false, // don't use steps (use length)
                                           false, // don't add paths
                                           go_forward,
                                           go_backward,
                                           id2);  // our target node
                    graph.rebuild_indexes();
#ifdef debug_mapper
                    if (debug) cerr << "got graph " << graph.size() << " " << pb2json(graph.graph) << endl;
#endif
                    //graph.serialize_to_file("raw-" + hash_alignment(aln) + ".vg");
                }

                // we have to remember how much we've trimmed from the first node
                // so that we can translate it after the fact
                map<id_t, pair<int, int> > trimmings;
                vector<id_t> target_nodes;

                // TODO continue if the graph doesn't have both cut points
                if (insertion_between_mems || !graph.has_node(id(first_cut)) || !graph.has_node(id(second_cut))) {
                    // treat the bit as unalignable
#ifdef debug_mapper
                    if (debug) cerr << "graph does not contain both cut points!" << endl;
#endif
                } else {

                    // now trim the graph to fit by cutting the head/tail node(s)
                    bool align_rc = false;
                    if (is_rev(first_cut) && is_rev(second_cut)) {
                        pos_t tmp_cut = first_cut;
                        first_cut = reverse(second_cut, graph.get_node(id(second_cut))->sequence().size());
                        second_cut = reverse(tmp_cut, graph.get_node(id(tmp_cut))->sequence().size());
                        align_rc = true;
                    } else {
                        if (is_rev(first_cut)) {
                            reverse(first_cut, graph.get_node(id(first_cut))->sequence().size());
                        }
                        if (is_rev(second_cut)) {
                            first_cut = reverse(second_cut, graph.get_node(id(second_cut))->sequence().size());
                        }
                    }

                    //cerr << "first_cut after " << first_cut << endl;
                    //cerr << "second_cut after " << second_cut << endl;

                    if (id(first_cut) == id(second_cut)) {
                        if (offset(first_cut) == offset(second_cut)) {
                            bool begin_cut = !offset(first_cut);
                            bool end_cut = (offset(first_cut) == graph.get_node(id(first_cut))->sequence().size());
                            if (!begin_cut && !end_cut) {
                                //cerr << "cut has offset" << endl;
                                Node* left; Node* right; Node* trimmed;
                                Node* node = graph.get_node(id(first_cut));
                                graph.divide_node(node, offset(first_cut), left, right);
                                //cerr << pb2json(*left) << " | " << pb2json(*right) << endl;
                                // check soft clip status, which will change what part we keep
                                // keep the part that's relevant to the soft clip resolution
                                if (soft_clip_to_left) {
                                    //cerr << "soft clip to left" << endl;
                                    graph.destroy_node(right);
                                    graph.swap_node_id(left, id(first_cut));
                                    trimmed = left;
                                    trimmings[id(first_cut)] = make_pair(0, offset(first_cut));
                                } else {
                                    //cerr << "soft clip to right or other" << endl;
                                    graph.destroy_node(left);
                                    graph.swap_node_id(right, id(first_cut));
                                    trimmed = right;
                                    trimmings[id(first_cut)] = make_pair(offset(first_cut), 0);
                                }
                                if (trimmed->sequence().size()) {
                                    target_nodes.push_back(trimmed->id());
                                } else {
                                    // push back each connected node
                                    for (auto& edge : graph.edges_to(trimmed)) {
                                        target_nodes.push_back(edge->from());
                                    }
                                    for (auto& edge : graph.edges_from(trimmed)) {
                                        target_nodes.push_back(edge->to());
                                    }
                                }
                            } else {
                                // erase everything before this node
                                // do so by removing edges
                                // later we will decide which subgraphs to keep
                                // check soft clip status, which will change what part we keep
                                if (soft_clip_to_left) {
                                    NodeSide keep = NodeSide(id(first_cut), false);
                                    for (auto& side : graph.sides_to(keep)) {
                                        target_nodes.push_back(side.node);
                                        graph.destroy_edge(side, keep);
                                    }
                                } else if (soft_clip_to_right) {
                                    NodeSide keep = NodeSide(id(first_cut), true);
                                    for (auto& side : graph.sides_from(keep)) {
                                        target_nodes.push_back(side.node);
                                        graph.destroy_edge(keep, side);
                                    }
                                } else {
                                    if (begin_cut) {
                                        assert(false);
                                    }
                                    if (end_cut) {
                                        assert(false);
                                    }
                                }
                            }
                        } else {
                            //cerr << "offsets different same node" << endl;
                            vector<int> positions = { (int)offset(first_cut), (int)offset(second_cut) };
                            vector<Node*> parts;
                            Node* node = graph.get_node(id(first_cut));
                            size_t orig_len = node->sequence().size();
                            graph.divide_node(node, positions, parts);
                            // now remove the end parts
                            graph.destroy_node(parts.front());
                            graph.destroy_node(parts.back());
                            graph.swap_node_id(parts.at(1), id(first_cut));
                            target_nodes.push_back(id(first_cut));
                            trimmings[id(first_cut)] = make_pair(offset(first_cut),
                                                                 orig_len - offset(second_cut));
                        }
                    } else { // different nodes to trim
                        //cerr << "different nodes" << endl;
                        if (offset(first_cut)) {
                            Node* left; Node* right;
                            Node* node = graph.get_node(id(first_cut));
                            graph.divide_node(node, offset(first_cut), left, right);
                            // remove the unused part
                            graph.destroy_node(left);
                            graph.swap_node_id(right, id(first_cut));
                            //target_nodes.push_back(graph.get_node(id(first_cut)));
                            Node* trimmed = graph.get_node(id(first_cut));
                            trimmings[id(first_cut)] = make_pair(offset(first_cut),
                                                                 0);
                            if (trimmed->sequence().size()) {
                                target_nodes.push_back(trimmed->id());
                            } else {
                                // push back each connected node
                                for (auto& edge : graph.edges_to(trimmed)) {
                                    target_nodes.push_back(edge->from());
                                }
                                for (auto& edge : graph.edges_from(trimmed)) {
                                    target_nodes.push_back(edge->to());
                                }
                            }
                        } else {
                            // destroy everything ahead of the node
                            NodeSide begin = NodeSide(id(first_cut));
                            for (auto& side : graph.sides_to(begin)) {
                                graph.destroy_edge(side, begin);
                            }
                            target_nodes.push_back(id(first_cut));
                        }

                        if (offset(second_cut)) {
                            Node* left; Node* right;
                            Node* node = graph.get_node(id(second_cut));
                            graph.divide_node(node, offset(second_cut), left, right);
                            // remove the unused part
                            graph.destroy_node(right);
                            graph.swap_node_id(left, id(second_cut));
                            //target_nodes.push_back(graph.get_node(id(second_cut)));
                            Node* trimmed = graph.get_node(id(first_cut));
                            if (trimmed->sequence().size()) {
                                target_nodes.push_back(trimmed->id());
                            } else {
                                // push back each connected node
                                for (auto& edge : graph.edges_to(trimmed)) {
                                    target_nodes.push_back(edge->from());
                                }
                                for (auto& edge : graph.edges_from(trimmed)) {
                                    target_nodes.push_back(edge->to());
                                }
                            }
                        } else {
                            // but we need to record the things in the graph connected to it
                            for (auto& side : graph.sides_to(id(second_cut))) {
                                target_nodes.push_back(side.node);
                            }
                            // destroy the node
                            graph.destroy_node(id(second_cut));
                            // we don't record this node as a target as we've destroyed it
                        }
                    }
                    graph.remove_null_nodes_forwarding_edges();
                    graph.remove_orphan_edges();
                }
                // reselect the target subgraph
                VG target;
                for (auto& id : target_nodes) {
#ifdef debug_mapper
                    //if (debug) cerr << "target node " << id << endl;
#endif
                    if (graph.has_node(id)) {
                        target.add_node(*graph.get_node(id));
                    }
                }
                graph.expand_context(target, edit.sequence().size(), false);
                graph = target;
                // now do the alignment
                if (graph.empty()) {
#ifdef debug_mapper
                    if (debug) {
                        cerr << "no target for alignment of " << edit.sequence()
                             << ", graph is empty" << endl;
                    }
#endif
                    score -= aligner->gap_open + edit.to_length()*aligner->gap_extension;
                    *new_mapping->add_edit() = edit;
                } else {
                    // we've set the graph to the trimmed target
#ifdef debug_mapper
                    if (debug) cerr << "target graph " << graph.size() << " " << pb2json(graph.graph) << endl;
#endif
                    //time to try an alignment
                    Alignment patch;
                    bool flip = mapping.position().is_reverse();
                    if (flip) {
                        patch.set_sequence(reverse_complement(edit.sequence()));
                        if (!aln.quality().empty()) {
                            string qual = aln.quality().substr(read_pos, edit.to_length());
                            reverse(qual.begin(), qual.end());
                            patch.set_quality(qual);
                        }
                    } else {
                        patch.set_sequence(edit.sequence());
                        if (!aln.quality().empty()) {
                            patch.set_quality(aln.quality().substr(read_pos, edit.to_length()));
                        }
                    }

                    // do the alignment
                    bool banded_global = !soft_clip_to_right && !soft_clip_to_left;
                    bool pinned_alignment = soft_clip_to_right || soft_clip_to_left;
                    bool pinned_reverse = false;
                    if (soft_clip_to_right) {
                        pinned_reverse = true;
                    }

                    patch = align_to_graph(patch,
                                           graph,
                                           max_query_graph_ratio,
                                           pinned_alignment,
                                           pinned_reverse,
                                           full_length_alignment_bonus,
                                           banded_global);

                    // adjust the translated node positions
                    for (int k = 0; k < patch.path().mapping_size(); ++k) {
                        auto* mapping = patch.mutable_path()->mutable_mapping(k);
                        auto t = trimmings.find(mapping->position().node_id());
                        if (t != trimmings.end()) {
                            auto trimmed_length_fwd = t->second.first;
                            auto trimmed_length_rev = t->second.second;
                            mapping->mutable_position()->set_offset(
                                mapping->position().offset() +
                                ( mapping->position().is_reverse() ? trimmed_length_rev : trimmed_length_fwd ));
                        }
                    }

                    // reverse complement back if we've flipped the read for alignment
                    if (flip) {
                        patch = reverse_complement_alignment(patch,
                                                             (function<int64_t(int64_t)>) ([&](int64_t id) {
                                                                     return (int64_t)get_node_length(id);
                                                                 }));
                    }

                    if (debug && !check_alignment(patch)) {
                        cerr << "patching failure " << pb2json(patched) << endl;
                        assert(false);
                    }

                    // append the chunk to patched

#ifdef debug_mapper
                    if (debug) cerr << "patch: " << pb2json(patch) << endl;
#endif
                    patch.clear_sequence(); // we set the whole sequence later
                    if (min_identity && patch.identity() < min_identity) {
                        //cerr << "doing that other thing" << endl;
                        score -= aligner->gap_open + edit.to_length()*aligner->gap_extension;
                        *new_mapping->add_edit() = edit;
                    } else {
                        //cerr << "extending alignment" << endl;
                        auto last_mapping = patched.mutable_path()->mutable_mapping(patched.path().mapping_size()-1);
                        if (last_mapping->edit_size() == 0
                            && last_mapping->position().node_id() != 0) {
                            // if we just did an alignment, use its position rather than a previous hint
                            // such as for soft clips
                            patched = merge_alignments(patch, patched, debug);
                        } else {
                            extend_alignment(patched, patch, true);

                        }
                        // point at the correct "new mapping"
                        new_mapping = patched.mutable_path()->mutable_mapping(patched.path().mapping_size()-1);
                        score += patch.score();
                    }
                    //cerr << "extended " << pb2json(patched) << endl;
                }
            }
            // update our offsets
            get_offset(ref_pos) += edit.from_length();
            read_pos += edit.to_length();
        }
        //cerr << "growing patched: " << pb2json(patched) << endl;
        /*
        #ifdef debug_mapper
    if (debug) {
            patched.set_sequence(aln.sequence().substr(0, read_pos));
            if (!check_alignment(patched)) {
                cerr << "patched failure " << pb2json(patched) << endl;
                assert(false);
            }
        }
        */
    }
    // finally, fix up the alignment score
    patched.set_sequence(aln.sequence());
    if (!aln.quality().empty()) {
        patched.set_quality(aln.quality());
    }

    // enable these to attempt alignment normalization
    patched = smooth_alignment(simplify(patched));
    patched.set_identity(identity(patched.path()));
    patched.set_score(score_alignment(patched));
    return patched;
}

Alignment Mapper::smooth_alignment(const Alignment& aln) {
    // find cases where we have reversals
    auto& path = aln.path();
    Alignment head, tail, smoothed;
    bool should_smooth = false;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        size_t to_len = mapping_to_length(mapping);
        size_t from_len = mapping_to_length(mapping);
        if (mapping_to_length(mapping) != mapping_from_length(mapping)) {
            should_smooth = true;
            break;
        }
        // two mappings to the same node
        if (i < path.mapping_size()-1
            && path.mapping(i).position().node_id()
            == path.mapping(i+1).position().node_id()) {
            should_smooth = true;
            break;
        }
    }
    if (should_smooth) {
#ifdef debug_mapper
        if (debug) cerr << "smoothing" << endl;
#endif
        // get the subgraph overlapping the alignment
        VG graph;
        id_t id1 = aln.path().mapping(0).position().node_id();
        id_t id2 = aln.path().mapping(aln.path().mapping_size()-1).position().node_id();
        id_t id_from = id1;
        id_t id_to = id2;
        if (id1 > id2) {
            id_from = id2;
            id_to = id1;
        }
        bool go_forward = true;
        bool go_backward = false;
        xindex->get_id_range(id_from, id_from, graph.graph);
        xindex->expand_context(graph.graph,
                               aln.sequence().size(),
                               false, // don't use steps (use length)
                               false, // don't add paths
                               go_forward,
                               go_backward,
                               id_to);  // our target node
        graph.rebuild_indexes();
        // re-do the alignment
        // against the graph
        // always use the banded global mode
        smoothed = aln;
        bool flip = aln.path().mapping(0).position().is_reverse();
        if (flip) {
            smoothed.set_sequence(reverse_complement(aln.sequence()));
            if (!aln.quality().empty()) {
                string qual = aln.quality();
                reverse(qual.begin(), qual.end());
                smoothed.set_quality(qual);
            }
        }
        bool banded_global = false;
        id_t pinned_id = 0;
        bool pinned_reverse = false;
        smoothed = align_to_graph(smoothed,
                                  graph,
                                  max_query_graph_ratio,
                                  pinned_id,
                                  pinned_reverse,
                                  banded_global);
        if (flip) {
            smoothed = reverse_complement_alignment(smoothed,
                                                    (function<int64_t(int64_t)>) ([&](int64_t id) {
                                                            return (int64_t)get_node_length(id);
                                                        }));
        }
        return simplify(smoothed);
    } else {
        return aln;
    }
}

// generate a score from the alignment without realigning
// handles split alignments, where gaps of unknown length are
// by estimating length using the positional paths embedded in the graph
int32_t Mapper::score_alignment(const Alignment& aln) {
    int score = 0;
    int read_offset = 0;
    auto& path = aln.path();
    auto aligner = aln.quality().empty() ? get_regular_aligner() : get_qual_adj_aligner();
    auto qual_adj_aligner = (QualAdjAligner*) aligner;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        //cerr << "looking at mapping " << pb2json(mapping) << endl;
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            //cerr << "looking at edit " << pb2json(edit) << endl;
            if (edit_is_match(edit)) {
                if (!aln.quality().empty()) {
                    score += qual_adj_aligner->score_exact_match(
                        aln.sequence().substr(read_offset, edit.to_length()),
                        aln.quality().substr(read_offset, edit.to_length()));
                } else {
                    score += edit.from_length()*aligner->match;
                }
            } else if (edit_is_deletion(edit)) {
                score -= aligner->gap_open + edit.from_length()*aligner->gap_extension;
            } else if (edit_is_insertion(edit)) {
                // todo how do we score this qual adjusted?
                score -= aligner->gap_open + edit.to_length()*aligner->gap_extension;
            }
            read_offset += edit.to_length();
        }
        // score any intervening gaps in mappings using approximate distances
        if (i+1 < path.mapping_size()) {
            // what is the distance between the last position of this mapping
            // and the first of the next
            Position last_pos = mapping.position();
            last_pos.set_offset(last_pos.offset() + mapping_from_length(mapping));
            Position next_pos = path.mapping(i+1).position();
#ifdef debug_mapper
            if (debug) cerr << "gap: " << make_pos_t(last_pos) << " to " << make_pos_t(next_pos) << endl;
#endif
            int dist = graph_distance(make_pos_t(last_pos), make_pos_t(next_pos), aln.sequence().size());
            if (dist == aln.sequence().size()) {
#ifdef debug_mapper
                if (debug) cerr << "could not find distance to next target" << endl;
#endif
                if (xindex->path_count) {
                    // use the embedded paths
                    dist = abs(xindex->min_distance_in_paths(
                                   last_pos.node_id(), last_pos.is_reverse(), last_pos.offset(),
                                   next_pos.node_id(), next_pos.is_reverse(), next_pos.offset()));
#ifdef debug_mapper
                    if (debug) cerr << "found distance using embedded paths" << endl;
#endif
                } else {
                    // we failed to find a distance
                    if (i+2 < path.mapping_size()) {
                        next_pos = path.mapping(i+2).position();
                        dist = graph_distance(make_pos_t(last_pos), make_pos_t(next_pos), aln.sequence().size());
                        dist -= mapping_from_length(path.mapping(i+1));
                    }
#ifdef debug_mapper
                    if (debug) cerr << "found distance using next mapping" << endl;
#endif
                }
            }
            if (dist) dist -= 1;
            // if the distance is the max...
#ifdef debug_mapper
            if (debug) cerr << "distance from " << pb2json(last_pos) << " to " << pb2json(next_pos) << " is " << dist << endl;
#endif
            if (dist > 0) {
                score -= aligner->gap_open + dist * aligner->gap_extension;
            }
        }
    }
    if (!softclip_start(aln)) {
        score += full_length_alignment_bonus;
    }
    if (!softclip_end(aln)) {
        score += full_length_alignment_bonus;
    }
#ifdef debug_mapper
    if (debug) cerr << "score from score_alignment " << score << endl;
#endif
    return max(0, score);
}

// make a perfect-match alignment out of a vector of MEMs which each have only one recorded hit
// use the base alignment sequence (which the SMEMs relate to) to fill in the gaps
Alignment Mapper::mems_to_alignment(const Alignment& aln, vector<MaximalExactMatch>& mems) {
    // base case--- empty alignment
    if (mems.empty()) {
        Alignment aln; return aln;
    }
    vector<Alignment> alns;
    // get reference to the start and end of the sequences
    string::const_iterator seq_begin = aln.sequence().begin();
    string::const_iterator seq_end = aln.sequence().end();
    // we use this to track where we need to add sequence
    string::const_iterator last_end = seq_begin;
    for (int i = 0; i < mems.size(); ++i) {
        auto& mem = mems.at(i);
        //cerr << "looking at " << mem.sequence() << endl;
        // this mem is contained in the last
        if (mem.end <= last_end) {
            continue;
        }
        // handle unaligned portion between here and the last SMEM or start of read
        if (mem.begin > last_end) {
            alns.emplace_back();
            alns.back().set_sequence(aln.sequence().substr(last_end - seq_begin, mem.begin - last_end));
        }
        Alignment aln = mem_to_alignment(mem);
        // find and trim overlap with previous
        if (i > 0) {
            // use the end of the last mem we touched (we may have skipped several)
            int overlap = last_end - mem.begin;
            if (overlap > 0) {
                aln = strip_from_start(aln, overlap);
            }
        }
        alns.push_back(aln);
        last_end = mem.end;
    }
    // handle unaligned portion at end of read
    alns.emplace_back();
    int start = last_end - seq_begin;
    int length = seq_end - (seq_begin + start);
    alns.back().set_sequence(aln.sequence().substr(start, length));
    auto alnm = merge_alignments(alns);
    *alnm.mutable_quality() = aln.quality();
    return alnm;
}

// convert one mem into an alignment; validates that only one node is given
Alignment Mapper::mem_to_alignment(MaximalExactMatch& mem) {
    const string seq = mem.sequence();
    if (mem.nodes.size() > 1) {
        cerr << "[vg::Mapper] warning: generating first alignment from MEM with multiple recorded hits" << endl;
    }
    auto& node = mem.nodes.front();
    pos_t pos = make_pos_t(gcsa::Node::id(node),
                           gcsa::Node::rc(node),
                           gcsa::Node::offset(node));
    return walk_match(seq, pos);
}

bool Mapper::get_mem_hits_if_under_max(MaximalExactMatch& mem) {
    bool filled = false;
    // remove all-Ns
    //if (!allATGC(mem.sequence())) return false;
    // require a minimum length
    if (mem.end-mem.begin == 0
        || mem.end-mem.begin < min_mem_length) return false;
    // use the counting interface to determine the number of hits
    mem.fill_match_count(gcsa);
    // if we aren't filtering on hit count, or if we have up to the max allowed hits
    if (mem.match_count > 0 && (!hit_max || mem.match_count <= hit_max)) {
        // extract the graph positions matching the range
        mem.fill_nodes(gcsa);
        filled = true;
    }
    return filled;
}

vector<Alignment> Mapper::align_mem_multi(const Alignment& alignment, vector<MaximalExactMatch>& mems, double& cluster_mq, int additional_multimaps) {

#ifdef debug_mapper
    if (debug) cerr << "aligning " << pb2json(alignment) << endl;
#endif
    if (!gcsa || !xindex) {
        cerr << "error:[vg::Mapper] a GCSA2/xg index pair is required for MEM mapping" << endl;
        exit(1);
    }

    if (mem_threading) {
        return mems_pos_clusters_to_alignments(alignment, mems, additional_multimaps, cluster_mq);
    } else {
        return mems_id_clusters_to_alignments(alignment, mems, additional_multimaps);
    }

}

vector<Alignment>
Mapper::mems_id_clusters_to_alignments(const Alignment& alignment, vector<MaximalExactMatch>& mems, int additional_multimaps) {

    struct StrandCounts {
        uint32_t forward;
        uint32_t reverse;
    };
    
    int total_multimaps = max_multimaps + additional_multimaps;

    // we will use these to determine the alignment strand for each subgraph
    map<id_t, StrandCounts> node_strands;
    // records a mapping of id->MEMs, for cluster ranking
    map<id_t, vector<MaximalExactMatch*> > id_to_mems;
    // for clustering
    vector<id_t> ids;

    // run through the mems, generating a set of alignments for each
    for (auto& mem : mems) {
        //#ifdef debug_mapper
        //if (debug) cerr << "on mem " << mem.sequence() << endl;
        //#endif
        size_t len = mem.begin - mem.end;
        // collect ids and orientations of hits to them on the forward mem
        for (auto& node : mem.nodes) {
            id_t id = gcsa::Node::id(node);
            id_to_mems[id].push_back(&mem);
            ids.push_back(id);
            if (gcsa::Node::rc(node)) {
                node_strands[id].reverse++;
            } else {
                node_strands[id].forward++;
            }
        }
    }

    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

    // establish clusters using approximate distance metric based on ids
    // we pick up ranges between successive nodes
    // when these are below our thread_extension length
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

    // rank the clusters by the fraction of the read that they cover
    map<vector<id_t>*, int> cluster_query_coverage;
    std::for_each(clusters.begin(), clusters.end(),
                  [&cluster_query_coverage,
                   &id_to_mems](vector<id_t>& cluster) {
                      set<string::const_iterator> query_coverage;
                      for (auto& id : cluster) {
                          auto& mems = id_to_mems[id];
                          std::for_each(mems.begin(), mems.end(),
                                        [&](MaximalExactMatch* m) {
                                            string::const_iterator c = m->begin;
                                            while (c != m->end) query_coverage.insert(c++);
                                        });
                      }
                      cluster_query_coverage[&cluster] = query_coverage.size();
                  });

    vector<vector<id_t>*> ranked_clusters;
    std::for_each(clusters.begin(), clusters.end(),
                  [&ranked_clusters](vector<id_t>& cluster) {
                      ranked_clusters.push_back(&cluster); });

    std::sort(ranked_clusters.begin(), ranked_clusters.end(),
              [&cluster_query_coverage](vector<id_t>* a,
                                        vector<id_t>* b) {
                  auto len_a = cluster_query_coverage[a];
                  auto len_b = cluster_query_coverage[b];
                  // order by cluster coverage of query
                  // break ties on number of MEMs (fewer better)
                  if (len_a == len_b) {
                      return a->size() < b->size();
                  } else {
                      return len_a > len_b;
                  }
              });


    // generate an alignment for each subgraph/orientation combination for which we have hits
#ifdef debug_mapper
    if (debug) cerr << "aligning to " << clusters.size() << " clusters" << endl;
    if (debug) {
        for (auto cptr : ranked_clusters) {
            auto& c = *cptr;
            cerr << cluster_query_coverage[cptr] << ":"
                 << c.size() << " "
                 << c.front() << "-" << c.back() << endl;
        }
    }
#endif

    vector<Alignment> alns; // our alignments
    
    // set up our forward and reverse base alignments (these are just sequences in bare alignment objs)
    auto aln_fw = alignment;
    aln_fw.clear_path();
    aln_fw.set_score(0);
    auto aln_rc = reverse_complement_alignment(aln_fw, (function<int64_t(int64_t)>)
                                               ([&](int64_t id) { return get_node_length(id); }));
    
    int max_target_length = alignment.sequence().size() * max_target_factor;

    size_t attempts = 0;
    for (auto& cptr : ranked_clusters) {
        auto& cluster = *cptr;
        // skip if our cluster is too small
        if (cluster.size() < cluster_min) continue;
        // record our attempt count
        ++attempts;
        // bail out if we've passed our maximum number of attempts
        if (attempts > max(max_attempts, total_multimaps)) break;
#ifdef debug_mapper
        if (debug) {
            cerr << "attempt " << attempts
                 << " on cluster " << cluster.front() << "-" << cluster.back() << endl;
        }
#endif
        VG sub; // the subgraph we'll align against
        set<id_t> seen;
        for (auto& id : cluster) {
            if (seen.count(id)) continue; // avoid double-gets
            seen.insert(id);
            xindex->get_id_range(id, id, sub.graph);
        }
        // expand using our context depth
        xindex->expand_context(sub.graph, context_depth, false);
        sub.rebuild_indexes();
        // if the graph is now too big to attempt, bail out
        if (max_target_factor && sub.length() > max_target_length) continue;
#ifdef debug_mapper
        if (debug) {
            cerr << "attempt " << attempts
                 << " on subgraph " << sub.min_node_id() << "-" << sub.max_node_id() << endl;
        }
#endif
        // determine the likely orientation
        uint32_t fw_mems = 0;
        uint32_t rc_mems = 0;
        sub.for_each_node([&](Node* n) {
                auto ns = node_strands.find(n->id());
                if (ns != node_strands.end()) {
                    fw_mems += ns->second.forward;
                    rc_mems += ns->second.reverse;
                }
            });
#ifdef debug_mapper
        if (debug) cerr << "got " << fw_mems << " forward and " << rc_mems << " reverse mems" << endl;
#endif
        if (fw_mems) {
            Alignment aln = align_to_graph(aln_fw, sub, max_query_graph_ratio);
            resolve_softclips(aln, sub);
            alns.push_back(aln);
            if (attempts >= total_multimaps &&
                greedy_accept &&
                aln.identity() >= accept_identity) {
                break;
            }
        }
        if (rc_mems) {
            Alignment aln = align_to_graph(aln_rc, sub, max_query_graph_ratio);
            resolve_softclips(aln, sub);
            alns.push_back(reverse_complement_alignment(aln,
                                                        (function<int64_t(int64_t)>)
                                                        ([&](int64_t id) { return get_node_length(id); })));
            if (attempts >= total_multimaps &&
                greedy_accept &&
                aln.identity() >= accept_identity) {
                break;
            }
        }
    }
    
    return alns;
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
    int max_target_length = aln.sequence().size() * max_target_factor;
    while (itr++ < max_softclip_iterations
           && (sc_start > softclip_threshold
               || sc_end > softclip_threshold)) {
#ifdef debug_mapper
        if (debug) {
            cerr << "Softclip before expansion: " << sc_start << " " << sc_end
                << " (" << aln.score() << " points)" << endl;
        }
#endif
        double avg_node_size = graph.length() / (double)graph.size();
#ifdef debug_mapper
        if (debug) cerr << "average node size " << avg_node_size << endl;
#endif
        // step towards the side where there were soft clips
        Graph flanks;
        xindex->get_id_range(idf, idf, flanks);
        xindex->get_id_range(idl, idl, flanks);
        xindex->expand_context(flanks,
                               max(context_depth, (int)((sc_start+sc_end)/avg_node_size)),
                               true, // use steps
                               false); // don't add paths
        graph.extend(flanks);

        aln.clear_path();
        aln.set_score(0);

        // give up if the graph is too big
        if (max_target_factor && graph.length() >= max_target_length) break;

        // otherwise, align
        aln = align_to_graph(aln, graph, max_query_graph_ratio);

        sc_start = softclip_start(aln);
        sc_end = softclip_end(aln);
#ifdef debug_mapper
        if (debug) {
            cerr << "Softclip after expansion: " << sc_start << " " << sc_end
                << " (" << aln.score() << " points)" << endl;
        }
#endif
        // we are not improving, so increasing the window is unlikely to help
        if (last_score == aln.score()) break;
        // update tracking of path end
        last_score = aln.score();
        path = aln.mutable_path();
        idf = path->mutable_mapping(0)->position().node_id();
        idl = path->mutable_mapping(path->mapping_size()-1)->position().node_id();
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
        //#ifdef debug_mapper
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
#ifdef debug_mapper
        if (debug) cerr << k << "\t~" << approx_matches << endl;
#endif
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
#ifdef debug_mapper
        if (debug) cerr << "\t=" << kmer_positions.size() << endl;
#endif
        kmer_count += kmer_positions.size();
        ++i;
    }

#ifdef debug_mapper
    if (debug) cerr << "kept kmer hits " << kmer_count << endl;
#endif

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

#ifdef debug_mapper
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
#endif

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
    #ifdef debug_mapper
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
    #endif
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

#ifdef debug_mapper
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
#endif

    int thread_ex = thread_extension;
    map<vector<int64_t>*, Alignment> alignments;
    auto aligner = (alignment.quality().empty() ? get_regular_aligner() : get_qual_adj_aligner());

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
#ifdef debug_mapper
            if (debug) cerr << "getting node range " << first << "-" << last << endl;
#endif
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

#ifdef debug_mapper
            if (debug) cerr << "got subgraph with " << graph->node_count() << " nodes, "
                            << graph->edge_count() << " edges" << endl;
#endif
            //serialize_to_file("init-" + alignment.sequence() + "-" + hash_alignment(alignment).substr(0,8) + "-" + hash().substr(0,8) + ".vg");

            // Topologically sort the graph, breaking cycles and orienting all edges end to start.
            // This flips some nodes around, so we need to translate alignments back.
            //set<int64_t> flipped_nodes;
            //graph->orient_nodes_forward(flipped_nodes);

            // align
            //graph->serialize_to_file("align2.vg");
            ta.clear_path();
            ta.set_score(0);

            ta = align_to_graph(ta, *graph, max_query_graph_ratio);

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
#ifdef debug_mapper
                if (debug) {
                    cerr << "softclip before " << sc_start << " " << sc_end << endl;
                    cerr << "distance to head "
                         << graph->distance_to_head(NodeTraversal(graph->get_node(idf), false), sc_start*3)
                         << endl;
                    cerr << "distance to tail "
                         << graph->distance_to_tail(NodeTraversal(graph->get_node(idl), false), sc_end*3)
                         << endl;
                }
#endif
                double avg_node_size = graph->length() / (double)graph->size();
#ifdef debug_mapper
                if (debug) cerr << "average node size " << avg_node_size << endl;
#endif
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

                ta = align_to_graph(ta, *graph, max_query_graph_ratio);

                sc_start = softclip_start(ta);
                sc_end = softclip_end(ta);
#ifdef debug_mapper
                if (debug) cerr << "softclip after " << sc_start << " " << sc_end << endl;
#endif
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

#ifdef debug_mapper
            if (debug) cerr << "normalized score is " << (float)ta.score() / ((float)ta.sequence().size()*(aligner->match)) << endl;
#endif
            if (greedy_accept && ta.identity() >= accept_identity) {
#ifdef debug_mapper
                if (debug) cerr << "greedy accept" << endl;
#endif
                accepted = true;
                break;
            }
        }
    }
    
    vector<Alignment> alns;
    alns.reserve(alignments.size());
    for (auto& ta : alignments) {
        alns.push_back(ta.second);
    }
    return alns;

}


// transform the path into a path relative to another path (defined by path_name)
// source -> surjection (in path_name coordinate space)
// the product is equivalent to a pairwise alignment between this path and the other

// new approach
// get path sequence
// get graph component overlapping path
// removing elements which aren't in the path of interest
// realign to this graph
// cross fingers

Alignment Mapper::surject_alignment(const Alignment& source,
                                    set<string>& path_names,
                                    string& path_name,
                                    int64_t& path_pos,
                                    bool& path_reverse,
                                    int window) {

    Alignment surjection = source;
    surjection.clear_mapping_quality();
    surjection.clear_score();
    surjection.clear_identity();
    surjection.clear_path();

    // get start and end nodes in path
    // get range between +/- window
    if (!source.has_path() || source.path().mapping_size() == 0) {
#ifdef debug

#pragma omp critical (cerr)
        cerr << "Alignment " << source.name() << " is unmapped and cannot be surjected" << endl;

#endif
        return surjection;
    }

    set<id_t> nodes;
    for (int i = 0; i < source.path().mapping_size(); ++ i) {
        nodes.insert(source.path().mapping(i).position().node_id());
    }
    VG graph;
    for (auto& node : nodes) {
        *graph.graph.add_node() = xindex->node(node);
    }
    xindex->expand_context(graph.graph, context_depth, true); // get connected edges and path
    graph.paths.append(graph.graph);
    graph.rebuild_indexes();

    set<string> kept_paths;
    graph.keep_paths(path_names, kept_paths);

    // We need this for inverting mappings to the correct strand
    function<int64_t(id_t)> node_length = [&graph](id_t node) {
        return graph.get_node(node)->sequence().size();
    };
    
    // What is our alignment to surject spelled the other way around? We can't
    // just use the normal alignment RC function because the mappings reference
    // nonexistent nodes.
    // Make sure to copy all the things about the alignment (name, etc.)

    Alignment surjection_rc = surjection;
    surjection_rc.set_sequence(reverse_complement(surjection.sequence()));
    
    // Align the old alignment to the graph in both orientations. Apparently
    // align only does a single oriantation, and we have no idea, even looking
    // at the mappings, which of the orientations will correspond to the one the
    // alignment is actually in.

    auto surjection_forward = align_to_graph(surjection, graph, max_query_graph_ratio);
    auto surjection_reverse = align_to_graph(surjection_rc, graph, max_query_graph_ratio);

#ifdef debug
#pragma omp critical (cerr)
    cerr << surjection.name() << " " << surjection_forward.score() << " forward score, " << surjection_reverse.score() << " reverse score" << endl;
#endif
    
    if(surjection_reverse.score() > surjection_forward.score()) {
        // Even if we have to surject backwards, we have to send the same string out as we got in.
        surjection = reverse_complement_alignment(surjection_reverse, node_length);
    } else {
        surjection = surjection_forward;
    }
    
    
#ifdef debug

#pragma omp critical (cerr)
        cerr << surjection.path().mapping_size() << " mappings, " << kept_paths.size() << " paths" << endl;

#endif

    if (surjection.path().mapping_size() > 0 && kept_paths.size() == 1) {
        // determine the paths of the node we mapped into
        //  ... get the id of the first node, get the paths of it
        assert(kept_paths.size() == 1);
        path_name = *kept_paths.begin();

        int64_t path_id = xindex->path_rank(path_name);
        auto& first_pos = surjection.path().mapping(0).position();
        int64_t hit_id = surjection.path().mapping(0).position().node_id();
        bool hit_backward = surjection.path().mapping(0).position().is_reverse();
        // we pick up positional information using the index

        auto path_posns = xindex->position_in_path(hit_id, path_name);
        if (path_posns.size() > 1) {
            cerr << "[vg map] surject_alignment: warning, multiple positions for node " << hit_id << " in " << path_name << " but will use only first: " << path_posns.front() << endl;
        } else if (path_posns.size() == 0) {
            cerr << "[vg map] surject_alignment: error, no positions for alignment " << source.name() << endl;
            exit(1);
        }

        // if we are reversed
        path_pos = path_posns.front();
        bool reversed_path = xindex->mapping_at_path_position(path_name, path_pos).position().is_reverse();
        if (reversed_path) {
            // if we got the start of the node position relative to the path
            // we need to offset to make thinsg right
            // but which direction
            if (hit_backward) {
                path_pos = path_posns.front() + first_pos.offset();
            } else {
                auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
            }
            path_reverse = !hit_backward;
        } else {
            if (!hit_backward) {
                path_pos = path_posns.front() + first_pos.offset();
            } else {
                auto pos = reverse_complement_alignment(surjection, node_length).path().mapping(0).position();
                path_pos = xindex->position_in_path(pos.node_id(), path_name).front() + pos.offset();
            }
            path_reverse = hit_backward;
        }

    } else {

        surjection = source;
#ifdef debug

#pragma omp critical (cerr)
        cerr << "Alignment " << source.name() << " did not align to the surjection subgraph" << endl;

#endif

    }

    return surjection;
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

bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin == m2.begin && m1.end == m2.end && m1.nodes == m2.nodes;
}

bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2) {
    return m1.begin < m2.begin && m1.end < m2.end && m1.nodes < m2.nodes;
}

MEMMarkovModel::MEMMarkovModel(
    const vector<size_t>& aln_lengths,
    const vector<vector<MaximalExactMatch> >& matches,
    Mapper* mapper,
    const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
    int band_width) {
    // store the MEMs in the model
    int frag_n = 0;
    for (auto& fragment : matches) {
        ++frag_n;
        for (auto& mem : fragment) {
            // copy the MEM for each specific hit in the base graph
            // and add it in as a vertex
            for (auto& node : mem.nodes) {
                //model.emplace_back();
                //auto m = model.back();
                MEMMarkovModelVertex m;
                m.weight = (double) mem.length() / (double)aln_lengths[frag_n-1];
                m.prev = nullptr;
                m.score = 0;
                m.mem = mem;
                m.mem.nodes.clear();
                m.mem.nodes.push_back(node);
                m.mem.fragment = frag_n;
                //m.mem.fill_positions(mapper);
                model.push_back(m);
            }
        }
    }
    for (vector<MEMMarkovModelVertex>::iterator m = model.begin(); m != model.end(); ++m) {
        // fill the nexts using banding constraints
        // banding means we stop looking along the MEMs when the score_transition function returns value/false
        // as a result changing the score function to return 0 at a given threshold induces local banding
        auto n = m;
        ++n;
        int i = 0;
        bool connected = false;
        //cerr << "from " << m->mem << endl;
        while (n != model.end()) {
            //cerr << "  to " << n->mem << endl;
            // skip past MEMs at the same position in the read
            if (n->mem.begin == m->mem.begin) {
                //cerr << "    skip" << endl;
                ++n; continue;
                if (n == model.end()) break;
            }
            // todo how do we handle MEMs at the same starting position
            // these are duplicates which we need to introduce...
            double weight = transition_weight(m->mem, n->mem);
            //cerr << "    weight " << weight << endl;
            if (++i > band_width && connected) {
                //cerr << "    breaking" << endl;
                break;
            }
            if (weight >= 0) {
                //cerr << "    saving" << endl;
                // save if we got a weight
                m->next_cost.push_back(make_pair(&*n, weight));
                n->prev_cost.push_back(make_pair(&*m, weight));
                connected = true;
            }
            ++n;
            ++i;
        }
        // sort the nexts to make later traversal easier
        sort(m->next_cost.begin(), m->next_cost.end(),
             [](const pair<MEMMarkovModelVertex*, double>& x,
                const pair<MEMMarkovModelVertex*, double>& y)
             { return x.second < y.second; });
    }
    //set<MEMMarkovModelVertex*> exclude; // empty set
    //score(exclude);
}

void MEMMarkovModel::score(const set<MEMMarkovModelVertex*>& exclude) {
    // propagate the scores in the model
    for (auto& m : model) {
        // score is equal to the max inbound + mem.weight
        if (exclude.count(&m)) continue; // skip
        m.score = m.weight; // base case
        for (auto& p : m.prev_cost) {
            double proposal = m.weight + p.second + p.first->score;
            if (proposal > m.score) {
                m.prev = p.first;
                m.score = proposal;
            }
        }
    }
}

MEMMarkovModelVertex* MEMMarkovModel::max_vertex(void) {
    MEMMarkovModelVertex* maxv = nullptr;
    for (auto& m : model) {
        if (maxv == nullptr || m.score > maxv->score) {
            maxv = &m;
        }
    }
    return maxv;
}

void MEMMarkovModel::clear_scores(void) {
    for (auto& m : model) {
        m.score = 0;
        m.prev = nullptr;
    }
}

vector<vector<MaximalExactMatch> > MEMMarkovModel::traceback(int alt_alns, bool debug) {
    set<MEMMarkovModelVertex*> used;
    vector<vector<MaximalExactMatch> > traces;
    for (int i = 0; i < alt_alns; ++i) {
        // score the model, accounting for excluded traces
        clear_scores();
        score(used);
#ifdef debug_mapper
        if (debug) {
            cerr << "MEMMarkovModel::traceback " << i << endl;
            display(cerr);
        }
#endif
        // find the maximum score
        auto* vertex = max_vertex();
        // check if we've exhausted our MEMs
        if (vertex == nullptr || vertex->score == 0) break;
        //cerr << "is maximum " << vertex->mem.sequence() << " " << vertex << ":" << vertex->score << endl;
        // make trace
        vector<MEMMarkovModelVertex*> vertex_trace;
        while (vertex != nullptr) {
            used.insert(vertex);
            vertex_trace.push_back(vertex);
            if (vertex->prev != nullptr) {
                vertex = vertex->prev;
            } else {
                break;
            }
        }
        traces.emplace_back();
        auto& mem_trace = traces.back();
        for (auto v = vertex_trace.rbegin(); v != vertex_trace.rend(); ++v) {
            mem_trace.push_back((*v)->mem);
        }
    }
    return traces;
}

// show model
void MEMMarkovModel::display(ostream& out) {
    for (auto& vertex : model) {
        out << vertex.mem.sequence() << ":" << vertex.mem.fragment << " " << &vertex << ":" << vertex.score << "@";
        for (auto& node : vertex.mem.nodes) {
            id_t id = gcsa::Node::id(node);
            size_t offset = gcsa::Node::offset(node);
            bool is_rev = gcsa::Node::rc(node);
            out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
        }
        out << "prev: ";
        for (auto& p : vertex.prev_cost) {
            auto& next = p.first;
            out << p.first << ":" << p.second << "@";
            for (auto& node : next->mem.nodes) {
                id_t id = gcsa::Node::id(node);
                size_t offset = gcsa::Node::offset(node);
                bool is_rev = gcsa::Node::rc(node);
                out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
            }
            out << " ; ";
        }
        out << " next: ";
        for (auto& p : vertex.next_cost) {
            auto& next = p.first;
            out << p.first << ":" << p.second << "@";
            for (auto& node : next->mem.nodes) {
                id_t id = gcsa::Node::id(node);
                size_t offset = gcsa::Node::offset(node);
                bool is_rev = gcsa::Node::rc(node);
                out << id << (is_rev ? "-" : "+") << ":" << offset << " ";
            }
            out << " ; ";
        }
        out << endl;
    }
}

// construct the sequence of the MEM; useful in debugging
string MaximalExactMatch::sequence(void) const {
    string seq; //seq.resize(end-begin);
    string::const_iterator c = begin;
    while (c != end) seq += *c++;
    return seq;
}

// uses GCSA to get the positions matching the range
void MaximalExactMatch::fill_nodes(gcsa::GCSA* gcsa) {
    gcsa->locate(range, nodes);
}

// uses GCSA to get a count of the number of graph nodes in our range
// if the match count is not already filled
void MaximalExactMatch::fill_match_count(gcsa::GCSA* gcsa) {
    if (match_count == 0) {
        match_count = gcsa->count(range);
    }
}

// length of the MEM
int MaximalExactMatch::length(void) const {
    return end - begin;
}

// uses an xgindex to fill out the MEM positions
void MaximalExactMatch::fill_positions(Mapper* mapper) {
    for (auto& node : nodes) {
        positions = mapper->node_positions_in_paths(gcsa::Node::encode(gcsa::Node::id(node), 0, gcsa::Node::rc(node)));
    }
}

ostream& operator<<(ostream& out, const MaximalExactMatch& mem) {
    size_t len = mem.begin - mem.end;
    out << mem.sequence() << ":";
    for (auto& node : mem.nodes) {
        id_t id = gcsa::Node::id(node);
        size_t offset = gcsa::Node::offset(node);
        bool is_rev = gcsa::Node::rc(node);
        out << id << (is_rev ? "-" : "+") << ":" << offset << ",";
    }
    return out;
}

}
