#include "alignment_io.hpp"
#include "gafkluge.hpp"

#include <sstream>
#include <regex>

namespace vg {

size_t unpaired_for_each_parallel(function<bool(Alignment&)> get_read_if_available,
                                  function<void(Alignment&)> lambda,
                                  uint64_t batch_size) {
    assert(batch_size % 2 == 0);    
    size_t nLines = 0;
    vector<Alignment> *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_read_if_available, lambda, batch_size)
#pragma omp single
    {
        
        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = batch_size;
        // max # we will ever increase the batch buffer to
        const uint64_t max_max_batches_outstanding = 1 << 13; // 8192
        
        // alignments to hold the incoming data
        Alignment aln;
        // did we find the end of the file yet?
        bool more_data = true;
        
        while (more_data) {
            // init a new batch
            batch = new std::vector<Alignment>();
            batch->reserve(batch_size);
            
            // load up to the batch-size number of reads
            for (int i = 0; i < batch_size; i++) {
                
                more_data = get_read_if_available(aln);
                
                if (more_data) {
                    batch->emplace_back(std::move(aln));
                    nLines++;
                }
                else {
                    break;
                }
            }
            
            // did we get a batch?
            if (batch->size()) {
                
                // how many batch tasks are outstanding currently, including this one?
                uint64_t current_batches_outstanding;
#pragma omp atomic capture
                current_batches_outstanding = ++batches_outstanding;
                
                if (current_batches_outstanding >= max_batches_outstanding) {
                    // do this batch in the current thread because we've spawned the maximum number of
                    // concurrent batch tasks
                    for (auto& aln : *batch) {
                        lambda(aln);
                    }
                    delete batch;
#pragma omp atomic capture
                    current_batches_outstanding = --batches_outstanding;
                    
                    if (4 * current_batches_outstanding / 3 < max_batches_outstanding
                        && max_batches_outstanding < max_max_batches_outstanding) {
                        // we went through at least 1/4 of the batch buffer while we were doing this thread's batch
                        // this looks risky, since we want the batch buffer to stay populated the entire time we're
                        // occupying this thread on compute, so let's increase the batch buffer size
                        
                        max_batches_outstanding *= 2;
                    }
                }
                else {
                    // spawn a new task to take care of this batch
#pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda)
                    {
                        for (auto& aln : *batch) {
                            lambda(aln);
                        }
                        delete batch;
#pragma omp atomic update
                        batches_outstanding--;
                    }
                }
            }
        }
    }
    return nLines;
}

size_t paired_for_each_parallel_after_wait(function<bool(Alignment&, Alignment&)> get_pair_if_available,
                                           function<void(Alignment&, Alignment&)> lambda,
                                           function<bool(void)> single_threaded_until_true,
                                           uint64_t batch_size) {

    assert(batch_size % 2 == 0);
    size_t nLines = 0;
    vector<pair<Alignment, Alignment> > *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
    
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_pair_if_available, single_threaded_until_true, lambda, batch_size)
#pragma omp single
    {

        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = batch_size;
        // max # we will ever increase the batch buffer to
        const uint64_t max_max_batches_outstanding = 1 << 13; // 8192
        
        // alignments to hold the incoming data
        Alignment mate1, mate2;
        // did we find the end of the file yet?
        bool more_data = true;
        
        while (more_data) {
            // init a new batch
            batch = new std::vector<pair<Alignment, Alignment>>();
            batch->reserve(batch_size);
            
            // load up to the batch-size number of pairs
            for (int i = 0; i < batch_size; i++) {
                
                more_data = get_pair_if_available(mate1, mate2);
                
                if (more_data) {
                    batch->emplace_back(std::move(mate1), std::move(mate2));
                    nLines++;
                }
                else {
                    break;
                }
            }
            
            // did we get a batch?
            if (batch->size()) {
                // how many batch tasks are outstanding currently, including this one?
                uint64_t current_batches_outstanding;
#pragma omp atomic capture
                current_batches_outstanding = ++batches_outstanding;
                
                bool do_single_threaded = !single_threaded_until_true();
                if (current_batches_outstanding >= max_batches_outstanding || do_single_threaded) {
                    // do this batch in the current thread because we've spawned the maximum number of
                    // concurrent batch tasks or because we are directed to work in a single thread
                    for (auto& p : *batch) {
                        lambda(p.first, p.second);
                    }
                    delete batch;
#pragma omp atomic capture
                    current_batches_outstanding = --batches_outstanding;
                    
                    if (4 * current_batches_outstanding / 3 < max_batches_outstanding
                        && max_batches_outstanding < max_max_batches_outstanding
                        && !do_single_threaded) {
                        // we went through at least 1/4 of the batch buffer while we were doing this thread's batch
                        // this looks risky, since we want the batch buffer to stay populated the entire time we're
                        // occupying this thread on compute, so let's increase the batch buffer size
                        // (skip this adjustment if you're in single-threaded mode and thus expect the buffer to be
                        // empty)
                        
                        max_batches_outstanding *= 2;
                    }
                }
                else {
                    // spawn a new task to take care of this batch
#pragma omp task default(none) firstprivate(batch) shared(batches_outstanding, lambda)
                    {
                        for (auto& p : *batch) {
                            lambda(p.first, p.second);
                        }
                        delete batch;
#pragma omp atomic update
                        batches_outstanding--;
                    }
                }
            }
        }
    }
    
    return nLines;
}

bool get_next_alignment_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer, gafkluge::GafRecord& g_buffer,
                                 Alignment& alignment) {
    if (hts_getline(fp, '\n', &s_buffer) <= 0) {
        return false;
    }

    gafkluge::parse_gaf_record(ks_str(&s_buffer), g_buffer);
    gaf_to_alignment(graph, g_buffer, alignment);
    return true;
}

bool get_next_interleaved_alignment_pair_from_gaf(const HandleGraph& graph, htsFile* fp, kstring_t& s_buffer,
                                                  gafkluge::GafRecord& g_buffer, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_gaf(graph, fp, s_buffer, g_buffer, mate1) &&
        get_next_alignment_from_gaf(graph, fp, s_buffer, g_buffer, mate2);
}

size_t gaf_unpaired_for_each(const HandleGraph& graph, const string& filename, function<void(Alignment&)> lambda) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln;
    gafkluge::GafRecord gaf;
    size_t count = 0;

    while (get_next_alignment_from_gaf(graph, in, s_buffer, gaf, aln) == true) {
        lambda(aln);
        ++count;
    }
    
    hts_close(in);

    return count;
}

size_t gaf_paired_interleaved_for_each(const HandleGraph& graph, const string& filename,
                                       function<void(Alignment&, Alignment&)> lambda) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    size_t count = 0;

    while (get_next_interleaved_alignment_pair_from_gaf(graph, in, s_buffer, gaf, aln1, aln2) == true) {
        lambda(aln1, aln2);
        count += 2;
    }
    
    hts_close(in);

    return count;
}

size_t gaf_unpaired_for_each_parallel(const HandleGraph& graph, const string& filename,
                                      function<void(Alignment&)> lambda,
                                      uint64_t batch_size) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }

    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    
    function<bool(Alignment&)> get_read = [&](Alignment& aln) {
        return get_next_alignment_from_gaf(graph, in, s_buffer, gaf, aln);
    };
        
    size_t nLines = unpaired_for_each_parallel(get_read, lambda, batch_size);
    
    hts_close(in);
    return nLines;

}

size_t gaf_paired_interleaved_for_each_parallel(const HandleGraph& graph, const string& filename,
                                                function<void(Alignment&, Alignment&)> lambda,
                                                uint64_t batch_size) {
    return gaf_paired_interleaved_for_each_parallel_after_wait(graph, filename, lambda, [](void) {return true;}, batch_size);
}

size_t gaf_paired_interleaved_for_each_parallel_after_wait(const HandleGraph& graph, const string& filename,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true,
                                                           uint64_t batch_size) {
    
    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }

    kstring_t s_buffer = KS_INITIALIZE;
    Alignment aln1, aln2;
    gafkluge::GafRecord gaf;
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_interleaved_alignment_pair_from_gaf(graph, in, s_buffer, gaf, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true, batch_size);

    hts_close(in);
    return nLines;    
}

gafkluge::GafRecord alignment_to_gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar, bool base_quals) {

    gafkluge::GafRecord gaf;

    //1 string Query sequence name
    gaf.query_name = aln.name();

    //2 int Query sequence length
    gaf.query_length = aln.sequence().length();

    //12 int Mapping quality (0-255; 255 for missing)
    //Note: protobuf can't distinguish between 0 and missing so we just copy it through
    gaf.mapq = aln.mapping_quality();

    if (aln.has_path() && aln.path().mapping_size() > 0) {    
        //3 int Query start (0-based; closed)
        gaf.query_start = 0; //(aln.path().mapping_size() ? first_path_position(aln.path()).offset() : "*") << "\t"
        //4 int Query end (0-based; open)
        gaf.query_end = aln.sequence().length();
        //5 char Strand relative to the path: "+" or "-"
        gaf.strand = '+'; // always positive relative to the path
        //7 int Path length
        gaf.path_length = 0;
        //8 int Start position on the path (0-based)
        gaf.path_start = 0;
        //10 int Number of residue matches
        gaf.matches = 0;
        gaf.path.reserve(aln.path().mapping_size());
        string cs_cigar_str;
        size_t running_match_length = 0;
        size_t total_to_len = 0;
        size_t prev_offset;
        for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
            auto& mapping = aln.path().mapping(i);
            size_t offset = mapping.position().offset();
            string node_seq;
            handle_t handle = graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
            bool skip_step = false;
            if (i == 0) {
                // use path_start to store the offset of the first node
                gaf.path_start = offset;
            } else if (cs_cigar == true && offset > 0) {
                if (offset == prev_offset && mapping.position().node_id() == aln.path().mapping(i -1).position().node_id() &&
                    mapping.position().is_reverse() == aln.path().mapping(i -1).position().is_reverse()) {
                    // our mapping is redundant, we won't write a step for it
                    skip_step = true;
                } else {
                    // to support split-mappings we gobble up the beginnings
                    // of nodes using deletions since unlike GAM, we can only
                    // set the offset of the first node
                    if (node_seq.empty()) {
                        node_seq = graph.get_sequence(handle);
                    }
                    cs_cigar_str += "-" + node_seq.substr(0, offset);
                }
            }
            for (size_t j = 0; j < mapping.edit_size(); ++j) {
                auto& edit = mapping.edit(j);
                if (edit_is_match(edit)) {
                    gaf.matches += edit.from_length();
                }
                if (cs_cigar == true) {
                    // CS-cigar string
                    if (edit_is_match(edit)) {
                        // Merge up matches that span edits/mappings
                        running_match_length += edit.from_length();
                    } else {
                        if (running_match_length > 0) {
                            // Matches are : followed by the match length
                            cs_cigar_str += ":" + std::to_string(running_match_length);
                            running_match_length = 0;
                        }
                        if (edit_is_sub(edit)) {
                            if (node_seq.empty()) {
                                node_seq = graph.get_sequence(handle);
                            }
                            // Substitions expressed one base at a time, preceded by *
                            for (size_t k = 0; k < edit.from_length(); ++k) {
                                cs_cigar_str += "*" + node_seq.substr(offset + k, 1) + edit.sequence().substr(k, 1); 
                            }
                        } else if (edit_is_deletion(edit)) {
                            if (node_seq.empty()) {
                                node_seq = graph.get_sequence(handle);
                            }
                            // Deletion is - followed by deleted sequence
                            assert(offset + edit.from_length() <= node_seq.length());
                            cs_cigar_str += "-" + node_seq.substr(offset, edit.from_length());
                        } else if (edit_is_insertion(edit)) {
                            // Insertion is "+" followed by inserted sequence
                            cs_cigar_str += "+" + edit.sequence();
                        }
                    }
                }
                offset += edit.from_length();
                total_to_len += edit.to_length();
            }

            if (i < aln.path().mapping_size() - 1 && offset != graph.get_length(handle)) {
                if (mapping.position().node_id() != aln.path().mapping(i + 1).position().node_id() ||
                    mapping.position().is_reverse() != aln.path().mapping(i + 1).position().is_reverse()) {
                    // we are hopping off the middle of a node, need to gobble it up with a deletion
                    if (node_seq.empty()) {
                        node_seq = graph.get_sequence(handle);
                    }
                    if (running_match_length > 0) {
                        // Matches are : followed by the match length
                        cs_cigar_str += ":" + std::to_string(running_match_length);
                        running_match_length = 0;
                    }
                    cs_cigar_str += "-" + node_seq.substr(offset);
                } else {
                    // we have a duplicate node mapping.  vg map actually produces these sometimes
                    // where an insert gets its own mapping even though its from_length is 0
                    // the gaf cigar format assumes nodes are fully covered, so we squish it out.
                    skip_step = true;
                }
            }
            
            //6 string Path matching /([><][^\s><]+(:\d+-\d+)?)+|([^\s><]+)/
            if (!skip_step) {
                auto& position = mapping.position();
                gafkluge::GafStep step;
                step.name = std::to_string(position.node_id());
                step.is_stable = false;
                step.is_reverse = position.is_reverse();
                step.is_interval = false;
                gaf.path_length += graph.get_length(graph.get_handle(position.node_id()));
                if (i == 0) {
                    gaf.path_start = position.offset();
                }
                gaf.path.push_back(std::move(step));
            }
            
            if (i == aln.path().mapping_size()-1) {
                //9 int End position on the path (0-based)
                gaf.path_end = offset;
                assert(gaf.path_end >= 0);
            }

            prev_offset = offset;
        }
        if (cs_cigar && running_match_length > 0) {
            cs_cigar_str += ":" + std::to_string(running_match_length);
            running_match_length = 0;
        }

        // We can support gam alignments without sequences by inferring the sequence length from edits
        if (gaf.query_length == 0 && total_to_len > 0) {
            gaf.query_length = total_to_len;
            gaf.query_end = total_to_len;
        } 

        //11 int Alignment block length
        gaf.block_length = std::max(gaf.path_end - gaf.path_start, gaf.query_length);

        // optional cs-cigar string
        if (cs_cigar) {
            gaf.opt_fields["cs"] = make_pair("Z", std::move(cs_cigar_str));
        }

        // convert the identity into the dv divergence field
        // https://lh3.github.io/minimap2/minimap2.html#10
        if (aln.identity() > 0) {
            stringstream dv_str;
            dv_str << std::floor((1. - aln.identity()) * 10000. + 0.5) / 10000.;
            gaf.opt_fields["dv"] = make_pair("f", dv_str.str());
        }

        // convert the score into the AS field
        // https://lh3.github.io/minimap2/minimap2.html#10
        if (aln.score() > 0) {
            gaf.opt_fields["AS"] = make_pair("i", std::to_string(aln.score()));
        }

        // optional base qualities
        if (base_quals && !aln.quality().empty()) { 
            gaf.opt_fields["bq"] = make_pair("Z", string_quality_short_to_char(aln.quality()));
        }   
                
    }

    return gaf;
    
}

void gaf_to_alignment(const HandleGraph& graph, const gafkluge::GafRecord& gaf, Alignment& aln) {

    aln.Clear();

    if (!gafkluge::is_missing(gaf.query_name)) {
        aln.set_name(gaf.query_name);
    }

    for (size_t i = 0; i < gaf.path.size(); ++i) {
        const auto& gaf_step = gaf.path[i];
        // only support unstable gaf at this point
        assert(gaf_step.is_stable == false);
        assert(gaf_step.is_interval == false);
        Mapping* mapping = aln.mutable_path()->add_mapping();
        mapping->mutable_position()->set_node_id(std::stol(gaf_step.name));
        mapping->mutable_position()->set_is_reverse(gaf_step.is_reverse);
        if (i == 0) {
            mapping->mutable_position()->set_offset(gaf.path_start);
        }
        mapping->set_rank(i + 1);
    }

    if (gaf.mapq != 255) {
        // We let 255 be equivalent to 0, which isn't great
        aln.set_mapping_quality(gaf.mapq);
    }

    if (!gaf.path.empty()) {
        size_t cur_mapping = 0;
        int64_t cur_offset = gaf.path_start;
        handle_t cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                               aln.path().mapping(cur_mapping).position().is_reverse());
        size_t cur_len = graph.get_length(cur_handle);
        string& sequence = *aln.mutable_sequence();
        // Use the CS cigar string to add Edits into our Path, as well as set the sequence
        gafkluge::for_each_cs(gaf, [&] (const string& cs_cigar) {
                assert(cur_offset < cur_len || (cs_cigar[0] == '+' && cur_offset <= cur_len));

                if (cs_cigar[0] == ':') {
                    int64_t match_len = stol(cs_cigar.substr(1));
                    while (match_len > 0) {
                        int64_t current_match = std::min(match_len, (int64_t)graph.get_length(cur_handle) - cur_offset);
                        Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                        edit->set_from_length(current_match);
                        edit->set_to_length(current_match);
                        sequence += graph.get_sequence(cur_handle).substr(cur_offset, current_match);
                        match_len -= current_match;
                        cur_offset += current_match;
                        if (match_len > 0) {
                            assert(cur_mapping < aln.path().mapping_size() - 1);
                            ++cur_mapping;
                            cur_offset = 0;
                            cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                                          aln.path().mapping(cur_mapping).position().is_reverse());
                            cur_len = graph.get_length(cur_handle);
                        }
                    }
                } else if (cs_cigar[0] == '+') {
                    size_t tgt_mapping = cur_mapping;
                    // left-align insertions to try to be more consistent with vg
                    if (cur_offset == 0 && cur_mapping > 0 && (!aln.path().mapping(cur_mapping - 1).position().is_reverse()
                                                               || cur_mapping == aln.path().mapping_size())) {
                        --tgt_mapping;
                    }
                    Edit* edit = aln.mutable_path()->mutable_mapping(tgt_mapping)->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(cs_cigar.length() - 1);
                    edit->set_sequence(cs_cigar.substr(1));
                    sequence += edit->sequence();
                } else if (cs_cigar[0] == '-') {
                    string del = cs_cigar.substr(1);
                    assert(del.length() <= graph.get_length(cur_handle) - cur_offset);
                    assert(del == graph.get_sequence(cur_handle).substr(cur_offset, del.length()));
                    Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                    edit->set_to_length(0);
                    edit->set_from_length(del.length());
                    cur_offset += del.length();
                    // unlike matches, we don't allow deletions to span multiple nodes
                    assert(cur_offset <= graph.get_length(cur_handle));
                } else if (cs_cigar[0] == '*') {
                    assert(cs_cigar.length() == 3);
                    char from = cs_cigar[1];
                    char to = cs_cigar[2];
                    assert(graph.get_sequence(cur_handle)[cur_offset] == from);
                    Edit* edit = aln.mutable_path()->mutable_mapping(cur_mapping)->add_edit();
                    // todo: support multibase snps
                    edit->set_from_length(1);
                    edit->set_to_length(1);
                    edit->set_sequence(string(1, to));
                    sequence += edit->sequence();
                    ++cur_offset;
                }
            
                // advance to the next mapping if we've pushed the offset past the current node
                assert(cur_offset <= cur_len);
                if (cur_offset == cur_len) {
                    ++cur_mapping;
                    cur_offset = 0;
                    if (cur_mapping < aln.path().mapping_size()) {
                        cur_handle = graph.get_handle(aln.path().mapping(cur_mapping).position().node_id(),
                                                      aln.path().mapping(cur_mapping).position().is_reverse());
                        cur_len = graph.get_length(cur_handle);
                    }
                }
            });
    }

    for (auto opt_it : gaf.opt_fields) {
        if (opt_it.first == "dv") {
            // get the identity from the dv divergence field
            // https://lh3.github.io/minimap2/minimap2.html#10
            aln.set_identity(1. - std::stof(opt_it.second.second));
        } else if (opt_it.first == "AS") {
            // get the score from the AS field
            // https://lh3.github.io/minimap2/minimap2.html#10
            aln.set_score(std::stoi(opt_it.second.second));
        } else if (opt_it.first == "bq") {
            // get the quality from the bq field
            aln.set_quality(string_quality_char_to_short(opt_it.second.second));
        }
    }
}

short quality_char_to_short(char c) {
    return static_cast<short>(c) - 33;
}

char quality_short_to_char(short i) {
    return static_cast<char>(i + 33);
}

void alignment_quality_short_to_char(Alignment& alignment) {
    alignment.set_quality(string_quality_short_to_char(alignment.quality()));
}

string string_quality_short_to_char(const string& quality) {
    string buffer; buffer.resize(quality.size());
    for (int i = 0; i < quality.size(); ++i) {
        buffer[i] = quality_short_to_char(quality[i]);
    }
    return buffer;
}

void alignment_quality_char_to_short(Alignment& alignment) {
    alignment.set_quality(string_quality_char_to_short(alignment.quality()));
}

string string_quality_char_to_short(const string& quality) {
    string buffer; buffer.resize(quality.size());
    for (int i = 0; i < quality.size(); ++i) {
        buffer[i] = quality_char_to_short(quality[i]);
    }
    return buffer;
}


}
