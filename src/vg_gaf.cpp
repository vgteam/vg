/**
 * \file vg_gaf.cpp
 *
 * GAM <==> GAF conversion
 */

#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "vg_gaf.hpp"
#include <htslib/hts.h>

namespace vg {
using namespace std;

gafkluge::GafRecord aln2gaf(const HandleGraph& graph, const Alignment& aln, bool cs_cigar) {

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
        //9 int End position on the path (0-based)
        uint64_t end_position = 0;
        //10 int Number of residue matches
        gaf.matches = 0;
        gaf.path.reserve(aln.path().mapping_size());
        string cs_cigar_str;
        size_t running_match_length = 0;
        size_t total_to_len = 0;
        for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
            auto& mapping = aln.path().mapping(i);
            size_t offset = mapping.position().offset();
            string node_seq;
            handle_t handle = graph.get_handle(mapping.position().node_id(), mapping.position().is_reverse());
            if (i == 0) {
                // use path_start to store the offset of the first node
                gaf.path_start = offset;
            } else if (cs_cigar == true && offset > 0) {
                // to support split-mappings we gobble up the beginnings
                // of nodes using deletions since unlike GAM, we can only
                // set the offset of the first node
                if (node_seq.empty()) {
                    node_seq = graph.get_sequence(handle);
                }
                cs_cigar_str += "-" + node_seq.substr(0, offset);
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

            bool skip_step = false;
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
                uint64_t node_length = graph.get_length(graph.get_handle(position.node_id()));
                gaf.path_length += node_length;
                if (i == 0) {
                    gaf.path_start = position.offset();
                }
                if (i == aln.path().mapping_size()-1) {
                    gaf.path_end = node_length - position.offset() - mapping_from_length(aln.path().mapping(i));
                }
                gaf.path.push_back(std::move(step));
            }
        }
        if (cs_cigar && running_match_length > 0) {
            cs_cigar_str += ":" + std::to_string(running_match_length);
            running_match_length = 0;
        }

        // We can support gam alignments without sequences by inferring the sequence length from edits
        if (gaf.query_length == 0 && total_to_len > 0) {
            gaf.query_length = total_to_len;
            gaf.query_end = total_to_len;
        } else if (gaf.query_length > 0 && total_to_len > 0 && gaf.query_length != total_to_len) {
            std::cerr << "Warning [gam2gaf]: Sequence length (" << gaf.query_length 
                      << ") not match total edit to length (" << total_to_len << " for " <<pb2json(aln) << endl;
        }

        //11 int Alignment block length
        gaf.block_length = std::max(gaf.path_end - gaf.path_start, gaf.query_length);

        // optional cs-cigar string
        if (cs_cigar) {
            gaf.opt_fields["cs"] = make_pair("Z", std::move(cs_cigar_str));
        }
    }

    return gaf;
    
}

Alignment gaf2aln(const HandleGraph& graph, const gafkluge::GafRecord& gaf) {

    Alignment aln;

    aln.set_name(gaf.query_name);

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
                assert(cur_offset < cur_len);

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
    return aln;
}

int gaf_for_each(const string& filename,
                 function<void(Alignment&)> lambda,
                 const HandleGraph& graph) {

    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        return 0;
    }
    
    kstring_t k_buffer = KS_INITIALIZE;
    Alignment aln;
    gafkluge::GafRecord gaf;
    while (hts_getline(in, '\n', &k_buffer) > 0) {
        
        gafkluge::parse_gaf_record(ks_str(&k_buffer), gaf);
        aln = gaf2aln(graph, gaf);
        lambda(aln);
    }
    hts_close(in);
    return 1;
}

static int gaf_for_each_parallel_impl(const string& filename,
                                      function<void(Alignment&, Alignment&)> lambda2,
                                      function<void(Alignment&)> lambda1, 
                                      const HandleGraph& graph) {

    // largely copied from hts_for_each_parallel in alignment.cpp
    
    // todo: this is a simplistic approach (all threads locking to read) that could probably be improved on by having a
    // reader thread filling up a queue.  ideally this would be done in a general framework to allow reading gaf/gam with
    // the same methods
    
    htsFile* in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "hts open fail on " << filename << endl;
        return 0;
    }

    bool more_data = true;    
    #pragma omp parallel shared(in, more_data)
    {
        kstring_t k_buffer1 = KS_INITIALIZE;
        kstring_t k_buffer2 = KS_INITIALIZE;                
        Alignment aln1;
        Alignment aln2;
        gafkluge::GafRecord gaf;

        while (more_data) {            
            // We need to track our own read operation's success separate from
            // the global flag, or someone else encountering EOF will cause us
            // to drop our read on the floor.
            bool got_read1 = false;
            bool got_read2 = false;
#pragma omp critical (hts_input)
            {
                if (more_data) {
                    got_read1 = hts_getline(in, '\n', &k_buffer1) > 0;
                    more_data &= got_read1;
                }
                if (more_data) {
                    got_read2 = hts_getline(in, '\n', &k_buffer2) > 0;
                    more_data &= got_read2;
                }
            }

            // Now we're outside the critical section so we can only rely on our own variables.
            if (got_read1) {
                gafkluge::parse_gaf_record(ks_str(&k_buffer1), gaf);
                aln1 = gaf2aln(graph, gaf);

                if (got_read2) {
                    gafkluge::parse_gaf_record(ks_str(&k_buffer2), gaf);
                    aln2 = gaf2aln(graph, gaf);
                    lambda2(aln1, aln2);
                } else {
                    lambda1(aln1);
                }
            }
        }
    }
    return 1;
}

int gaf_for_each_parallel(const string& filename,
                          function<void(Alignment&)> lambda,
                          const HandleGraph& graph) {

    function<void(Alignment&, Alignment&)> lambda2 = [&](Alignment& aln1, Alignment& aln2) {
        lambda(aln1);
        lambda(aln2);
    };

    return gaf_for_each_parallel_impl(filename, lambda2, lambda, graph);    
}

int gaf_for_each_interleaved_pair_parallel(const string& filename,
                                           function<void(Alignment&, Alignment&)> lambda2,
                                           const HandleGraph& graph) {
    function<void(Alignment&)> err1 = [](Alignment&){
        throw runtime_error("gaf_for_each_interleaved_pair_parallel: expected input stream of interleaved pairs, but it had odd number of elements");
    };
    
    return gaf_for_each_parallel_impl(filename, lambda2, err1, graph);
}


}
