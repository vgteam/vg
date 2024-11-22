#include "alignment.hpp"
#include "vg/io/gafkluge.hpp"
#include "annotation.hpp"
#include <vg/io/stream.hpp>

#include <sstream>
#include <chrono>

using namespace vg::io;

namespace vg {

int hts_for_each(string& filename, function<void(Alignment&)> lambda, const PathPositionHandleGraph* graph) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);
    map<int, path_handle_t> tid_path_handle;
    parse_tid_path_handle_map(hdr, graph, tid_path_handle);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        Alignment a = bam_to_alignment(b, rg_sample, tid_path_handle, hdr, graph);
        lambda(a);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

int hts_for_each(string& filename, function<void(Alignment&)> lambda) {
    return hts_for_each(filename, lambda, nullptr);
}

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda,
                          const PathPositionHandleGraph* graph) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);
    map<int, path_handle_t> tid_path_handle;
    parse_tid_path_handle_map(hdr, graph, tid_path_handle);

    int thread_count = get_thread_count();
    vector<bam1_t*> bs; bs.resize(thread_count);
    for (auto& b : bs) {
        b = bam_init1();
    }

    bool more_data = true;
#pragma omp parallel shared(in, hdr, more_data, rg_sample)
    {
        int tid = omp_get_thread_num();
        while (more_data) {
            bam1_t* b = bs[tid];
            // We need to track our own read operation's success separate from
            // the global flag, or someone else encountering EOF will cause us
            // to drop our read on the floor.
            bool got_read = false;
#pragma omp critical (hts_input)
            if (more_data) {
                got_read = sam_read1(in, hdr, b) >= 0;
                more_data &= got_read;
            }
            // Now we're outside the critical section so we can only rely on our own variables.
            if (got_read) {
                Alignment a = bam_to_alignment(b, rg_sample, tid_path_handle, hdr, graph);
                lambda(a);
            }
        }
    }

    for (auto& b : bs) bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda) {
    return hts_for_each_parallel(filename, lambda, nullptr);
}

bam_hdr_t* hts_file_header(string& filename, string& header) {
    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) {
        cerr << "[vg::alignment] could not open " << filename << endl;
        exit(1);
    }
    bam_hdr_t *hdr = sam_hdr_read(in);
    header = hdr->text;
    bam_hdr_destroy(hdr);
    hts_close(in);
    return hdr;
}

bam_hdr_t* hts_string_header(string& header,
                             const map<string, int64_t>& path_length,
                             const map<string, string>& rg_sample) {
    
    // Copy the map into a vecotr in its own order
    vector<pair<string, int64_t>> path_order_and_length(path_length.begin(), path_length.end());
    
    // Make header in that order.
    return hts_string_header(header, path_order_and_length, rg_sample);
}

bam_hdr_t* hts_string_header(string& header,
                             const vector<pair<string, int64_t>>& path_order_and_length,
                             const map<string, string>& rg_sample) {
    stringstream hdr;
    hdr << "@HD\tVN:1.5\tSO:unknown\n";
    for (auto& p : path_order_and_length) {
        hdr << "@SQ\tSN:" << p.first << "\t" << "LN:" << p.second << "\n";
    }
    for (auto& s : rg_sample) {
        hdr << "@RG\tID:" << s.first << "\t" << "SM:" << s.second << "\n";
    }
    hdr << "@PG\tID:0\tPN:vg\n";
    header = hdr.str();
    string sam = "data:," + header;
    samFile *in = sam_open(sam.c_str(), "r");
    bam_hdr_t *h = sam_hdr_read(in);
    sam_close(in);
    return h;
}

bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment) {

    alignment.Clear();
    bool is_fasta = false;
    // handle name
    string name;
    if (gzgets(fp,buffer,len) != 0) {
        buffer[strlen(buffer)-1] = '\0';
        name = buffer;
        if (name[0] == '@') {
            is_fasta = false;
        } else if (name[0] == '>') {
            is_fasta = true;
        } else {
            throw runtime_error("Found unexpected delimiter " + name.substr(0,1) + " in fastq/fasta input");
        }
        name = name.substr(1, name.find(' ') - 1); // trim off leading @ and things after the first whitespace
        // keep trailing /1 /2
        alignment.set_name(name);
    }
    else {
        // no more to get
        return false;
    }
    // handle sequence
    string sequence;
    bool reading_sequence = true;
    while (reading_sequence) {
        if (gzgets(fp,buffer,len) == 0) {
            if (sequence.empty()) {
                // there was no sequence
                throw runtime_error("[vg::alignment.cpp] incomplete fastq/fasta record " + name);
            }
            else {
                // we hit the end of the file
                break;
            }
        }
        size_t size_read = strlen(buffer);
        if (buffer[size_read - 1] == '\n') {
            // we stopped because of a line end rather than because we filled the buffer
            
            // we don't want the newline in the sequence, so terminate the buffer 1 char earlier
            --size_read;
            if (!is_fasta) {
                // we assume FASTQ sequences only take one line
                reading_sequence = false;
            }
            else {
                // peek ahead to check for a multi-line sequence
                int c = gzgetc(fp);
                if (c < 0) {
                    // this is the end of the file
                    reading_sequence = false;
                }
                else {
                    if (c == '>') {
                        // the next line is a sequence name
                        reading_sequence = false;
                    }
                    // un-peek
                    gzungetc(c, fp);
                }
            }
        }
        sequence.append(buffer, size_read);
    }
    alignment.set_sequence(sequence);
    // handle "+" sep
    if (!is_fasta) {
        if (0!=gzgets(fp,buffer,len)) {
        } else {
            cerr << "[vg::alignment.cpp] error: incomplete fastq record " << name << endl; exit(1);
        }
        // handle quality
        if (0!=gzgets(fp,buffer,len)) {
            buffer[strlen(buffer)-1] = '\0';
            string quality = string_quality_char_to_short(buffer);
            //cerr << string_quality_short_to_char(quality) << endl;
            alignment.set_quality(quality);
        } else {
            cerr << "[vg::alignment.cpp] error: fastq record missing base quality " << name << endl; exit(1);
        }
    }

    return true;

}

bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp, buffer, len, mate1) && get_next_alignment_from_fastq(fp, buffer, len, mate2);
}

bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp1, buffer, len, mate1) && get_next_alignment_from_fastq(fp2, buffer, len, mate2);
}

size_t fastq_unpaired_for_each_parallel(const string& filename, function<void(Alignment&)> lambda, uint64_t batch_size) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 2 << 22; // 4M
    char* buf = new char[len];
    
    function<bool(Alignment&)> get_read = [&](Alignment& aln) {
        return get_next_alignment_from_fastq(fp, buf, len, aln);;
    };
    
    
    size_t nLines = unpaired_for_each_parallel(get_read, lambda, batch_size);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
    
}

size_t fastq_paired_interleaved_for_each_parallel(const string& filename, function<void(Alignment&, Alignment&)> lambda, uint64_t batch_size) {
    return fastq_paired_interleaved_for_each_parallel_after_wait(filename, lambda, [](void) {return true;}, batch_size);
}
    
size_t fastq_paired_two_files_for_each_parallel(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda, uint64_t batch_size) {
    return fastq_paired_two_files_for_each_parallel_after_wait(file1, file2, lambda, [](void) {return true;}, batch_size);
}
    
size_t fastq_paired_interleaved_for_each_parallel_after_wait(const string& filename,
                                                             function<void(Alignment&, Alignment&)> lambda,
                                                             function<bool(void)> single_threaded_until_true,
                                                             uint64_t batch_size) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 1 << 18; // 256k
    char* buf = new char[len];
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_interleaved_alignment_pair_from_fastq(fp, buf, len, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true, batch_size);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
}
    
size_t fastq_paired_two_files_for_each_parallel_after_wait(const string& file1, const string& file2,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true,
                                                           uint64_t batch_size) {
    
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp1) {
        cerr << "[vg::alignment.cpp] couldn't open " << file1 << endl; exit(1);
    }
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp2) {
        cerr << "[vg::alignment.cpp] couldn't open " << file2 << endl; exit(1);
    }
    
    size_t len = 1 << 18; // 256k
    char* buf = new char[len];
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_alignment_pair_from_fastqs(fp1, fp2, buf, len, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true, batch_size);
    
    delete[] buf;
    gzclose(fp1);
    gzclose(fp2);
    return nLines;
}

size_t fastq_unpaired_for_each(const string& filename, function<void(Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    size_t len = 2 << 22; // 4M
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment alignment;
    while(get_next_alignment_from_fastq(fp, buffer, len, alignment)) {
        lambda(alignment);
        nLines++;
    }
    gzclose(fp);
    delete[] buffer;
    return nLines;
}

size_t fastq_paired_interleaved_for_each(const string& filename, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment mate1, mate2;
    while(get_next_interleaved_alignment_pair_from_fastq(fp, buffer, len, mate1, mate2)) {
        lambda(mate1, mate2);
        nLines++;
    }
    gzclose(fp);
    delete[] buffer;
    return nLines;
}


size_t fastq_paired_two_files_for_each(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp1) {
        cerr << "[vg::alignment.cpp] couldn't open " << file1 << endl; exit(1);
    }
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp2) {
        cerr << "[vg::alignment.cpp] couldn't open " << file2 << endl; exit(1);
    }
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment mate1, mate2;
    while(get_next_alignment_pair_from_fastqs(fp1, fp2, buffer, len, mate1, mate2)) {
        lambda(mate1, mate2);
        nLines++;
    }
    gzclose(fp1);
    gzclose(fp2);
    delete[] buffer;
    return nLines;

}

void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample) {
    string header(hts_header);
    vector<string> header_lines = split_delims(header, "\n");

    for (auto& line : header_lines) {

        // get next line from header, skip if empty
        if ( line.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if (line.find("@RG") == 0) {
            vector<string> rg_parts = split_delims(line, "\t ");
            string name;
            string rg_id;
            for (auto& part : rg_parts) {
                size_t colpos = part.find(":");
                if (colpos != string::npos) {
                    string fieldname = part.substr(0, colpos);
                    if (fieldname == "SM") {
                        name = part.substr(colpos+1);
                    } else if (fieldname == "ID") {
                        rg_id = part.substr(colpos+1);
                    }
                }
            }
            if (name.empty()) {
                cerr << "[vg::alignment] Error: could not find 'SM' in @RG line " << endl << line << endl;
                exit(1);
            }
            if (rg_id.empty()) {
                cerr << "[vg::alignment] Error: could not find 'ID' in @RG line " << endl << line << endl;
                exit(1);
            }
            map<string, string>::iterator s = rg_sample.find(rg_id);
            if (s != rg_sample.end()) {
                if (s->second != name) {
                    cerr << "[vg::alignment] Error: multiple samples (SM) map to the same read group (RG)" << endl
                          << endl
                          << "samples " << name << " and " << s->second << " map to " << rg_id << endl
                          << endl
                          << "It will not be possible to determine what sample an alignment belongs to" << endl
                          << "at runtime." << endl
                          << endl
                          << "To resolve the issue, ensure that RG ids are unique to one sample" << endl
                          << "across all the input files to freebayes." << endl
                          << endl
                          << "See bamaddrg (https://github.com/ekg/bamaddrg) for a method which can" << endl
                          << "add RG tags to alignments." << endl;
                    exit(1);
                }
            }
            // if it's the same sample name and RG combo, no worries
            rg_sample[rg_id] = name;
        }
    }
}

void parse_tid_path_handle_map(const bam_hdr_t* hts_header, const PathHandleGraph* graph, map<int, path_handle_t>& tid_path_handle) {
    if (!graph) {
        // No path handles to find!
        return;
    }
    for (int i = 0; i < hts_header->n_targets; i++) {
        // Pre-look-up all the paths mentioned in the header
        string target_name(hts_header->target_name[i]);
        if (graph->has_path(target_name)) {
            path_handle_t target = graph->get_path_handle(target_name);
            if (graph->get_sense(target) != PathSense::HAPLOTYPE) {
                // Non-haplotype paths are allowed in the mapping because they
                // are always path-position indexed.
                
                // Store the handles for the paths we find, under their HTSlib target numbers.
                tid_path_handle.emplace(i, target);
            } else {
                // TODO: Decide we need to positional-index this path? Make
                // PackedReferencePathOverlay take a collection of paths to
                // index and use this one?
                #pragma omp critical (cerr)
                std::cerr << "error[vg::parse_tid_path_handle_map] Path " << target_name
                          << " referenced in header exists in graph, but as a haplotype."
                          << " It is probably not indexed for positional lookup. Make the"
                          << " path a reference path"
                          << " <https://github.com/vgteam/vg/wiki/Changing-References>"
                          << " and try again." << std::endl;
                exit(1);
            }
        }
    }
}

// Internal conversion function for both paired and unpaired codepaths
string alignment_to_sam_internal(const Alignment& alignment,
                                 const string& refseq,
                                 const int32_t refpos,
                                 const bool refrev,
                                 const vector<pair<int, char>>& cigar,
                                 const string& mateseq,
                                 const int32_t matepos,
                                 bool materev,
                                 const int32_t tlen,
                                 bool paired,
                                 const int32_t tlen_max) {

    // Determine flags, using orientation, next/prev fragments, and pairing status.
    int32_t flags = determine_flag(alignment, refseq, refpos, refrev, mateseq, matepos, materev, tlen, paired, tlen_max);
    
    string alignment_name;
    if (paired) {
        // We need to strip the /1 and /2 or _1 and _2 from paired reads so the two ends have the same name.
        alignment_name = regex_replace(alignment.name(), regex("[/_][12]$"), "");
    } else {
        // Keep the alignment name as is because even if the name looks paired, the reads are semantically unpaired.
        alignment_name = alignment.name();
    }
    
    // Have One True Flag for whether the read is mapped (and should have its
    // mapping stuff set) or unmapped (and should have things *'d out).
    bool mapped = !(flags & BAM_FUNMAP);
        
    if (mapped) {
        // Make sure we have everything
        assert(!refseq.empty());
        assert(refpos != -1);
        assert(!cigar.empty());
        assert(alignment.has_path());
        assert(alignment.path().mapping_size() > 0);
    }

    // We apply the convention of unmapped reads getting their mate's coordinates
    // See section 2.4.1 https://samtools.github.io/hts-specs/SAMv1.pdf
    bool use_mate_loc = !mapped && paired && !mateseq.empty();
    
    stringstream sam;
    
    sam << (!alignment_name.empty() ? alignment_name : "*") << "\t"
        << flags << "\t"
        << (mapped ? refseq : use_mate_loc ? mateseq : "*") << "\t"
        << (use_mate_loc ? matepos + 1 : refpos + 1) << "\t"
        << (mapped ? alignment.mapping_quality() : 0) << "\t"
        << (mapped ? cigar_string(cigar) : "*") << "\t"
        << (mateseq == "" ? "*" : (mateseq == refseq ? "=" : mateseq)) << "\t"
        << matepos + 1 << "\t"
        << tlen << "\t"
        // Make sure sequence always comes out in reference forward orientation by looking at the flags.
        << (!alignment.sequence().empty() ? (refrev ? reverse_complement(alignment.sequence()) : alignment.sequence()) : "*") << "\t";
    if (!alignment.quality().empty()) {
        auto quality = alignment.quality();
        if (refrev) {
            // Quality also needs to be flipped
            std::reverse(quality.begin(), quality.end());
        }
        for (int i = 0; i < quality.size(); ++i) {
            sam << quality_short_to_char(quality[i]);
        }
    } else {
        sam << "*";
    }
    //<< (alignment.has_quality() ? string_quality_short_to_char(alignment.quality()) : string(alignment.sequence().size(), 'I'));
    if (!alignment.read_group().empty()) sam << "\tRG:Z:" << alignment.read_group();
    if (has_annotation(alignment, "tags")) {
        for (const auto& tag : get_annotation<vector<string>>(alignment, "tags")) {
            sam << '\t' << tag;
        }
    }
    sam << "\n";
    return sam.str();
}

int32_t determine_flag(const Alignment& alignment,
                       const string& refseq,
                       const int32_t refpos,
                       const bool refrev,
                       const string& mateseq,
                       const int32_t matepos,
                       bool materev,
                       const int32_t tlen,
                       bool paired,
                       const int32_t tlen_max) {
    
    // Determine flags, using orientation, next/prev fragments, and pairing status.
    int32_t flags = sam_flag(alignment, refrev, paired);
    
    // We've observed some reads with the unmapped flag set and also a CIGAR string set, which shouldn't happen.
    // We will check for this. The CIGAR string will only be set in the output if the alignment has a path.
    assert((bool)(flags & BAM_FUNMAP) != (alignment.has_path() && alignment.path().mapping_size()));
    
    if (!((bool)(flags & BAM_FUNMAP)) && paired && !refseq.empty() && refseq == mateseq) {
        // Properly paired if both mates mapped to same sequence, in inward-facing orientations.
        // We know they're on the same sequence, so check orientation.
        
        // If we are first, mate needs to be reverse, and if mate is first, we need to be reverse.
        // If we are at the same position either way is fine.
        bool facing = ((refpos <= matepos) && !refrev && materev) || ((matepos <= refpos) && refrev && !materev);
        
        // We are close enough if there is not tlen limit, or if there is one and we do not exceed it
        bool close_enough = (tlen_max == 0) || abs(tlen) <= tlen_max;
        
        if (facing && close_enough) {
            // We can't find anything wrong with this pair; it's properly paired.
            flags |= BAM_FPROPER_PAIR;
        }
        
        // TODO: Support sequencing technologies where "proper" pairing may
        // have a different meaning or expected combination of orientations.
    }
    
    if (paired && mateseq.empty()) {
        // Set the flag for the mate being unmapped
        flags |= BAM_FMUNMAP;
    }
    
    if (paired && materev) {
        // Set the flag for the mate being reversed
        flags |= BAM_FMREVERSE;
    }
    
    return flags;
}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const vector<pair<int, char>>& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        bool materev,
                        const int32_t tlen,
                        const int32_t tlen_max) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, mateseq, matepos, materev, tlen, true, tlen_max);

}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const vector<pair<int, char>>& cigar) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, "", -1, false, 0, false, 0);

}

vector<tuple<string, char, string>> parse_sam_tags(const vector<string>& tags) {
    
    vector<tuple<string, char, string>> parsed;
    for (const auto& tag : tags) {
        if (tag.empty()) {
            continue;
        }
        if (tag.size() < 6 || tag[2] != ':' || tag[4] != ':') {
            std::cerr << ("error: failed to parse malformed SAM tag '" + tag + "'\n");
            exit(1);
        }
        parsed.emplace_back(tag.substr(0, 2), tag[3], tag.substr(5, string::npos));
    }
    return parsed;
}

// template to reduce redunant code parsing and writing B type SAM tags
template<typename T>
void write_array_to_aux(bam1_t* bam, const char* tag_name, const string& arr_string) {
    
    vector<T> parsed;
    for (const auto& token : split_delims(arr_string.substr(1, string::npos), ",")) {
        parsed.push_back(parse<T>(token));
    }
    // size includes array type and length
    size_t data_size = parsed.size() * sizeof(T) + 5;
    uint8_t* data = (uint8_t*) malloc(data_size);
    // add the type
    data[0] = arr_string[0];
    // add the length
    *((uint32_t*) (data + 1)) = (uint32_t) parsed.size();
    // add the array
    for (size_t i = 0, j = 5; i < parsed.size(); ++i, j += sizeof(T)) {
        *((T*) (data + j)) = parsed[i];
    }
    bam_aux_append(bam, tag_name, 'B', data_size, data);
    free(data);
}

// Internal conversion function for both paired and unpaired codepaths
bam1_t* alignment_to_bam_internal(bam_hdr_t* header,
                                  const Alignment& alignment,
                                  const string& refseq,
                                  const int32_t refpos,
                                  const bool refrev,
                                  const vector<pair<int, char>>& cigar,
                                  const string& mateseq,
                                  const int32_t matepos,
                                  bool materev,
                                  const int32_t tlen,
                                  bool paired,
                                  const int32_t tlen_max) {
    
    // this table doesn't seem to be reproduced in htslib publicly, so I'm copying
    // it from the CRAM conversion code
    static const char nt_encoding[256] = {
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15, 0,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15, 1,14, 2,13,15,15, 4,11,15,15,12,15, 3,15,15,
        15,15, 5, 6, 8,15, 7, 9,15,10,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
        15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15
    };
    
    // init an empty BAM record
    bam1_t* bam = bam_init1();
    
    // strip the pair order identifiers
    string alignment_name = alignment.name();
    if (paired && alignment_name.size() >= 2) {
        // We need to strip the /1 and /2 or _1 and _2 from paired reads so the two ends have the same name.
        char c1 = alignment_name[alignment_name.size() - 2];
        char c2 = alignment_name[alignment_name.size() - 1];
        if ((c1 == '_' || c1 == '/') && (c2 == '1' || c2 == '2')) {
            alignment_name = alignment_name.substr(0, alignment_name.size() - 2);
        }
    }
    
    // calculate the size in bytes of the variable length fields (which are all concatenated in memory)
    int qname_nulls = 4 - alignment_name.size() % 4;
    int qname_data_size = alignment_name.size() + qname_nulls;
    int cigar_data_size = 4 * cigar.size();
    int seq_data_size = (alignment.sequence().size() + 1) / 2; // round up
    int qual_data_size = alignment.sequence().size(); // we will allocate this even if quality doesn't exist
    
    // allocate the joint variable length fields
    int var_field_data_size = qname_data_size + cigar_data_size + seq_data_size + qual_data_size;
    bam->data = (uint8_t*) calloc(var_field_data_size, sizeof(uint8_t));
    
    // TODO: what ID is this? CRAM seems to ignore it, so maybe we can too...
    //bam->id = 0;
    bam->l_data = var_field_data_size; // current length of data
    bam->m_data = var_field_data_size; // max length of data
    
    bam1_core_t& core = bam->core;
    // mapping position
    core.pos = refpos;
    // ID of sequence mapped to
    core.tid = sam_hdr_name2tid(header, refseq.c_str());
    // MAPQ
    core.qual = alignment.mapping_quality();
    // number of nulls (above 1) used to pad read name string
    core.l_extranul = qname_nulls - 1;
    // bit flag
    core.flag = determine_flag(alignment, refseq, refpos, refrev, mateseq, matepos, materev, tlen, paired, tlen_max);
    // length of read name, including nulls
    core.l_qname = qname_data_size;
    // number of cigar operations
    core.n_cigar = cigar.size();
    // length of read
    core.l_qseq = alignment.sequence().size();
    // ID of sequence mate is mapped to
    core.mtid = sam_hdr_name2tid(header, mateseq.c_str()); // TODO: what if there is no mate
    // mapping position of mate
    core.mpos = matepos;
    // insert length of fragment
    core.isize = tlen;
    
    // all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    
    // write query name, padded by nulls
    uint8_t* name_data = bam->data;
    for (size_t i = 0; i < alignment_name.size(); ++i) {
        name_data[i] = (uint8_t) alignment_name[i];
    }
    for (size_t i = 0; i < qname_nulls; ++i) {
        name_data[i + alignment_name.size()] = '\0';
    }
    
    // encode cigar and copy into data

    uint32_t* cigar_data = (uint32_t*) (name_data + qname_data_size);
    
    auto refend = core.pos;
    for (size_t i = 0; i < cigar.size(); ++i) {
        uint32_t op;
        switch (cigar[i].second) {
            case 'M':
            case 'm':
                op = BAM_CMATCH;
                refend += cigar[i].first;
                break;
            case 'I':
            case 'i':
                op = BAM_CINS;
                break;
            case 'D':
            case 'd':
                op = BAM_CDEL;
                refend += cigar[i].first;
                break;
            case 'N':
            case 'n':
                op = BAM_CREF_SKIP;
                refend += cigar[i].first;
                break;
            case 'S':
            case 's':
                op = BAM_CSOFT_CLIP;
                break;
            case 'H':
            case 'h':
                op = BAM_CHARD_CLIP;
                break;
            case 'P':
            case 'p':
                op = BAM_CPAD;
                break;
            case '=':
                op = BAM_CEQUAL;
                refend += cigar[i].first;
                break;
            case 'X':
            case 'x':
                op = BAM_CDIFF;
                refend += cigar[i].first;
                break;
            default:
                throw runtime_error("Invalid CIGAR operation " + string(1, cigar[i].second));
                break;
        }
        cigar_data[i] = bam_cigar_gen(cigar[i].first, op);
    }
    
    
    // now we know where it ends, we can compute the bin
    // copied from cram/cram_samtools.h
    core.bin = hts_reg2bin(refpos, refend - 1, 14, 5); // TODO: not sure if end is past-the-last
    
    // convert sequence to 4-bit (nibble) encoding
    uint8_t* seq_data = (uint8_t*) (cigar_data + cigar.size());
    const string* seq = &alignment.sequence();
    string rev_seq;
    const string* qual = &alignment.quality();
    string rev_qual;
    if (refrev) {
        // Sequence and quality both need to be flipped to target forward orientation
        rev_seq = reverse_complement(*seq);
        seq = &rev_seq;
        reverse_copy(qual->begin(), qual->end(), back_inserter(rev_qual));
        qual = &rev_qual;
    }
    for (size_t i = 0; i < alignment.sequence().size(); i += 2) {
        if (i + 1 < alignment.sequence().size()) {
            seq_data[i / 2] = (nt_encoding[seq->at(i)] << 4) | nt_encoding[seq->at(i + 1)];
        }
        else {
            seq_data[i / 2] = nt_encoding[seq->at(i)] << 4;
        }
    }
    
    // write the quality directly (it should already have the +33 offset removed)
    uint8_t* qual_data = seq_data + seq_data_size;
    for (size_t i = 0; i < alignment.sequence().size(); ++i) {
        if (alignment.quality().empty()) {
            // hacky, but this seems to be what they do in CRAM anyway
            qual_data[i] = '\xff';
        }
        else {
            qual_data[i] = qual->at(i);
        }
    }
    
    if (~core.flag & BAM_FUNMAP) {
        // we've decided that it is aligned
        int32_t score = alignment.score();
        bam_aux_append(bam, "AS", 'i', sizeof(int32_t), (uint8_t*) &score);
    }
    
    if (!alignment.read_group().empty()) {
        bam_aux_append(bam, "RG", 'Z', alignment.read_group().size() + 1, (uint8_t*) alignment.read_group().c_str());
    }
    
    // this annotation comes from surject and should be retained in the BAM
    if (has_annotation(alignment, "all_scores")) {
        string all_scores = get_annotation<string>(alignment, "all_scores");
        bam_aux_append(bam, "SS", 'Z', all_scores.size() + 1, (uint8_t*) all_scores.c_str());
    }
    
    // TODO: it would be nice wrap htslib and set the other tags this way as well
    if (has_annotation(alignment, "tags")) {
        auto parsed_tags = parse_sam_tags(get_annotation<vector<string>>(alignment, "tags"));
        for (const auto& tag : parsed_tags) {
            
            const char* tag_id = get<0>(tag).c_str();
            char tag_type = get<1>(tag);
            const string& tag_val = get<2>(tag);
            if (get<0>(tag).size() != 2) {
                cerr << ("error: SAM tag label " + get<0>(tag) + " is not 2 characters long\n");
                exit(1);
            }
            if (tag_val.empty()) {
                cerr << ("error: SAM tag " + get<0>(tag) + " is missing a value\n");
                exit(1);
            }
            
            switch (tag_type) {
                case 'A':
                    // character
                    if (tag_val.size() != 1) {
                        cerr << ("error: SAM tag of type 'A' is not a single character: " + tag_val + "\n");
                        exit(1);
                    }
                    bam_aux_append(bam, tag_id, tag_type, sizeof(char), (uint8_t*) &tag_val[0]);
                    break;
                case 'i':
                    // integer
                {
                    int32_t val = parse<int32_t>(tag_val);
                    bam_aux_append(bam, tag_id, tag_type, sizeof(int32_t), (uint8_t*) &val);
                    break;
                }
                case 'f':
                    // float
                {
                    float val = parse<float>(tag_val);
                    bam_aux_append(bam, tag_id, tag_type, sizeof(float), (uint8_t*) &val);
                    break;
                }
                case 'Z':
                    // string
                case 'H':
                    // hex strings are copied as raw strings
                    bam_aux_append(bam, tag_id, tag_type, tag_val.size() + 1, (uint8_t*) tag_val.c_str());
                    break;
                case 'B':
                {
                    // the array of values has its own sub-type for entries
                    char subtype = tag_val.front();
                    switch (subtype) {
                        case 'c':
                            write_array_to_aux<int8_t>(bam, tag_id, tag_val);
                            break;
                        case 'C':
                            write_array_to_aux<uint8_t>(bam, tag_id, tag_val);
                            break;
                        case 's':
                            write_array_to_aux<int16_t>(bam, tag_id, tag_val);
                            break;
                        case 'S':
                            write_array_to_aux<uint16_t>(bam, tag_id, tag_val);
                            break;
                        case 'i':
                            write_array_to_aux<int32_t>(bam, tag_id, tag_val);
                            break;
                        case 'I':
                            write_array_to_aux<uint32_t>(bam, tag_id, tag_val);
                            break;
                        case 'f':
                            write_array_to_aux<float>(bam, tag_id, tag_val);
                            break;
                        default:
                            cerr << ("error: unrecognized array type '" + string(1, subtype) + "' in 'B' type SAM tag\n");
                            exit(1);
                            break;
                    }
                    break;
                }
                default:
                    cerr << ("error: unrecognized SAM tag type '" + string(1, tag_type) + "'\n");
                    exit(1);
                    break;
            }
        }
    }
    
    // TODO: this does not seem to be a standardized field (https://samtools.github.io/hts-specs/SAMtags.pdf)
//    if (!alignment.sample_name()) {
//
//    }
        
    return bam;
}

bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         bool materev,
                         const int32_t tlen,
                         const int32_t tlen_max) {

    return alignment_to_bam_internal(bam_header, alignment, refseq, refpos, refrev, cigar, mateseq, matepos, materev, tlen, true, tlen_max);

}

bam1_t* alignment_to_bam(bam_hdr_t* bam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const bool refrev,
                         const vector<pair<int, char>>& cigar) {
    
    return alignment_to_bam_internal(bam_header, alignment, refseq, refpos, refrev, cigar, "", -1, false, 0, false, 0);

}

string cigar_string(const vector<pair<int, char> >& cigar) {
    vector<pair<int, char> > cigar_comp;
    pair<int, char> cur = make_pair(0, '\0');
    for (auto& e : cigar) {
        if (cur == make_pair(0, '\0')) {
            cur = e;
        } else {
            if (cur.second == e.second) {
                cur.first += e.first;
            } else {
                cigar_comp.push_back(cur);
                cur = e;
            }
        }
    }
    cigar_comp.push_back(cur);
    stringstream cigarss;
    for (auto& e : cigar_comp) {
        cigarss << e.first << e.second;
    }
    return cigarss.str();
}

string mapping_string(const string& source, const Mapping& mapping) {
    string result;
    int p = mapping.position().offset();
    for (const auto& edit : mapping.edit()) {
        // mismatch/sub state
// *matches* from_length == to_length, or from_length > 0 and offset unset
// *snps* from_length == to_length; sequence = alt
        // mismatch/sub state
        if (edit.from_length() == edit.to_length()) {
            if (!edit.sequence().empty()) {
                result += edit.sequence();
            } else {
                result += source.substr(p, edit.from_length());
            }
            p += edit.from_length();
        } else if (edit.from_length() == 0 && edit.sequence().empty()) {
// *skip* from_length == 0, to_length > 0; implies "soft clip" or sequence skip
            //cigar.push_back(make_pair(edit.to_length(), 'S'));
        } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
            result += edit.sequence();
            p += edit.from_length();
        } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
            result += edit.sequence();
            p += edit.from_length();
        }
    }
    return result;
}

void mapping_cigar(const Mapping& mapping, vector<pair<int, char>>& cigar) {
    for (const auto& edit : mapping.edit()) {
        if (edit.from_length() && edit.from_length() == edit.to_length()) {
// *matches* from_length == to_length, or from_length > 0 and offset unset
            // match state
            append_cigar_operation(edit.from_length(), 'M', cigar);
            //cerr << "match " << edit.from_length() << endl;
        } else {
            // mismatch/sub state
// *snps* from_length == to_length; sequence = alt
            if (edit.from_length() == edit.to_length()) {
                append_cigar_operation(edit.from_length(), 'M', cigar);
                //cerr << "match " << edit.from_length() << endl;
            } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
                int32_t del = edit.from_length() - edit.to_length();
                int32_t eq = edit.to_length();
                if (eq) append_cigar_operation(eq, 'M', cigar);
                append_cigar_operation(del, 'D', cigar);
                //cerr << "del " << edit.from_length() - edit.to_length() << endl;
            } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
                int32_t ins = edit.to_length() - edit.from_length();
                int32_t eq = edit.from_length();
                if (eq) append_cigar_operation(eq, 'M', cigar);
                append_cigar_operation(ins, 'I', cigar);
                //cerr << "ins " << edit.to_length() - edit.from_length() << endl;
            }
        }
    }
}

int64_t cigar_mapping(const bam1_t *b, Mapping* mapping) {
    int64_t ref_length = 0;
    int64_t query_length = 0;

    const auto cigar = bam_get_cigar(b);

    for (int k = 0; k < b->core.n_cigar; k++) {
        Edit* e = mapping->add_edit();
        const int op = bam_cigar_op(cigar[k]);
        const int ol = bam_cigar_oplen(cigar[k]);
        if (bam_cigar_type(cigar[k])&1) {
            // Consume query
            e->set_to_length(ol);
            string sequence; sequence.resize(ol);
            for (int i = 0; i < ol; i++ ) {
               sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(b), query_length + i)];
            }
            e->set_sequence(sequence);
            query_length += ol;
        } else {
            e->set_to_length(0);
        }
        if (bam_cigar_type(cigar[k])&2) {
            // Consume ref
            e->set_from_length(ol);
            ref_length += ol;
        } else {
            e->set_from_length(0);
        }
    }
    return ref_length;
}

void mapping_against_path(Alignment& alignment, const bam1_t *b, const path_handle_t& path, const PathPositionHandleGraph* graph, bool on_reverse_strand) {

    if (b->core.pos == -1) return;

    Mapping mapping;

    int64_t length = cigar_mapping(b, &mapping);

    Alignment aln = target_alignment(graph, path, b->core.pos, b->core.pos + length, "", on_reverse_strand, mapping);

    *alignment.mutable_path() = aln.path();

    Position* refpos = alignment.add_refpos();
    refpos->set_name(graph->get_path_name(path));
    refpos->set_offset(b->core.pos);
    refpos->set_is_reverse(on_reverse_strand);
}

vector<pair<int, char>> cigar_against_path(const Alignment& alignment, bool on_reverse_strand, int64_t& pos, size_t path_len, size_t softclip_suppress) {
    vector<pair<int, char> > cigar;
    
    if (!alignment.has_path() || alignment.path().mapping_size() == 0) return cigar;
    const Path& path = alignment.path();
    int l = 0;

    for (const auto& mapping : path.mapping()) {
        mapping_cigar(mapping, cigar);
    }
    
    if(on_reverse_strand) {
        // Flip CIGAR ops into forward strand ordering
        reverse(cigar.begin(), cigar.end());
    }

    // handle soft clips, which are just insertions at the start or end
    // back
    if (cigar.size() > 1 && cigar.back().second == 'D' && cigar[cigar.size() - 2].second == 'I') {
        // Swap insert to the outside so it can be a softclip.
        // When making the CIGAR we should put D before I but when flipping the
        // strand they may switch.
        std::swap(cigar.back(), cigar[cigar.size() - 2]);
    }
    if (cigar.back().second == 'I') {
        // make sure we stay in the reference sequence when suppressing the softclips
        if (cigar.back().first <= softclip_suppress
            && pos + alignment_from_length(alignment) + cigar.back().first <= path_len) {
            cigar.back().second = 'M';
        } else {
            cigar.back().second = 'S';
        }
    }
    // front
    if (cigar.size() > 1 && cigar.front().second == 'D' && cigar[1].second == 'I') {
        // Swap insert to the outside so it can be a softclip
        std::swap(cigar.front(), cigar[1]);
    }
    if (cigar.front().second == 'I') {
        // make sure we stay in the reference sequence when suppressing the softclips
        if (cigar.front().first <= softclip_suppress
            && pos - cigar.front().first >= 0) {
            cigar.front().second = 'M';
            pos -= cigar.front().first;
        } else {
            cigar.front().second = 'S';
        }
    }
    
    simplify_cigar(cigar);

    return cigar;
}

void simplify_cigar(vector<pair<int, char>>& cigar) {
    
    size_t removed = 0;
    for (size_t i = 0, j = 0; i < cigar.size(); ++j) {
        if (j == cigar.size() || (cigar[j].second != 'I' && cigar[j].second != 'D')) {
            // this is the end boundary of a runs of I/D operations
            if (j - i >= 3) {
                // we have at least 3 adjacent I/D operations, which means they should
                // be re-consolidated
                int d_total = 0, i_total = 0;
                for (size_t k = i - removed, end = j - removed; k < end; ++k) {
                    if (cigar[k].second == 'D') {
                        d_total += cigar[k].first;
                    }
                    else {
                        i_total += cigar[k].first;
                    }
                }
                
                cigar[i - removed] = make_pair(d_total, 'D');
                cigar[i - removed + 1] = make_pair(i_total, 'I');
                
                // mark that we've removed cigar operations
                removed += j - i - 2;
            }
            // move the start of the next I/D run beyond the current operation
            i = j + 1;
        }
        if (j < cigar.size()) {
            cigar[j - removed] = cigar[j];
        }
    }
    cigar.resize(cigar.size() - removed);
    // do a second pass removing empty operations and consolidating non I/D operations
    removed = 0;
    for (size_t i = 0; i < cigar.size(); ++i) {
        if (cigar[i].first == 0) {
            ++removed;
        }
        else if (i > removed && cigar[i].second == cigar[i - removed - 1].second) {
            cigar[i - removed - 1].first += cigar[i].first;
            ++removed;
        }
        else if (removed) {
            cigar[i - removed] = cigar[i];
        }
    }
    cigar.resize(cigar.size() - removed);
}

pair<int32_t, int32_t> compute_template_lengths(const int64_t& pos1, const vector<pair<int, char>>& cigar1,
    const int64_t& pos2, const vector<pair<int, char>>& cigar2) {

    // Compute signed distance from outermost matched/mismatched base of each
    // alignment to the outermost matched/mismatched base of the other.
    
    // We work with CIGARs because it's easier than reverse complementing
    // Alignment objects without node lengths.
    
    // Work out the low and high mapped bases for each side
    auto find_bounds = [](const int64_t& pos, const vector<pair<int, char>>& cigar) {
        // Initialize bounds to represent no mapped bases
        int64_t low = numeric_limits<int64_t>::max();
        int64_t high = numeric_limits<int64_t>::min();
        
        // Track position in the reference
        int64_t here = pos;
        for (auto& item : cigar) {
            // Trace along the cigar
            if (item.second == 'M') {
                // Bases are matched. Count them in the bounds and execute the operation
                low = min(low, here);
                here += item.first;
                high = max(high, here);
            } else if (item.second == 'D') {
                // Only other way to advance in the reference
                here += item.first;
            }
        }
        
        return make_pair(low, high);
    };
    
    auto bounds1 = find_bounds(pos1, cigar1);
    auto bounds2 = find_bounds(pos2, cigar2);
    
    // Compute the separation
    int32_t dist = 0;
    if (bounds1.first < bounds2.second) {
        // The reads are in order
        dist = bounds2.second - bounds1.first;
    } else if (bounds2.first < bounds1.second) {
        // The reads are out of order so the other bounds apply
        dist = bounds1.second - bounds2.first;
    }
    
    if (pos1 < pos2) {
        // Count read 1 as the overall "leftmost", so its value will be positive
        return make_pair(dist, -dist);
    } else {
        // Count read 2 as the overall leftmost
        return make_pair(-dist, dist);
    }

}

int32_t sam_flag(const Alignment& alignment, bool on_reverse_strand, bool paired) {
    int16_t flag = 0;

    if (paired) {
        // Respect the alignment's internal crossreferences.
        // Allow for multiple-read-long fragments. 
        
        flag |= BAM_FPAIRED;
        if (!alignment.has_fragment_next()) {
            // This is the last read in a pair
            flag |= BAM_FREAD2;
        }
        if (!alignment.has_fragment_prev()) {
            // This is the first read in a pair
            flag |= BAM_FREAD1;
        }
        
        // Invalid paired GAM is caught, for surject, on GAM input
        // TODO: catch reads with pair partners when they shouldn't be paired?
    }

    if (!alignment.has_path() || alignment.path().mapping_size() == 0) {
        // unmapped
        flag |= BAM_FUNMAP;
    } 
    if (on_reverse_strand) {
        flag |= BAM_FREVERSE;
    }
    if (alignment.is_secondary()) {
        flag |= BAM_FSECONDARY;
    }
    
    return flag;
}

vector<string> bam_tag_strings(const bam1_t* b) {
    
    vector<string> tag_strings;
    for (uint8_t* iter = bam_aux_first(b); iter != NULL; iter = bam_aux_next(b, iter)) {
        // tag name
        char name[2];
        name[0] = iter[0];
        name[1] = iter[1];
        kstring_t kstr;
        ks_initialize(&kstr);
        bam_aux_get_str(b, name, &kstr);
        tag_strings.emplace_back(ks_c_str(&kstr));
        ks_release(&kstr);
    }
    return tag_strings;
}

Alignment bam_to_alignment(const bam1_t *b,
                           const map<string, string>& rg_sample,
                           const map<int, path_handle_t>& tid_path_handle,
                           const bam_hdr_t *bh,
                           const PathPositionHandleGraph* graph) {

    Alignment alignment;

    // get the sequence and qual
    int32_t lqseq = b->core.l_qseq;
    string sequence; sequence.resize(lqseq);

    uint8_t* qualptr = bam_get_qual(b);
    string quality;//(lqseq, 0);
    quality.assign((char*)qualptr, lqseq);

    // process the sequence into chars
    uint8_t* seqptr = bam_get_seq(b);
    for (int i = 0; i < lqseq; ++i) {
        sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }

    // Now name the read after the scaffold
    string read_name = bam_get_qname(b);

    // Decide if we are a first read (/1) or second (last) read (/2)
    if(b->core.flag & BAM_FREAD1) {
        read_name += "/1";
    }
    if(b->core.flag & BAM_FREAD2) {
        read_name += "/2";
    }
    
    // If we are marked as both first and last we get /1/2, and if we are marked
    // as neither the scaffold name comes through unchanged as the read name.
    // TODO: produce correct names for intermediate reads on >2 read scaffolds.

    // add features to the alignment
    alignment.set_name(read_name);
    // was the sequence reverse complemented?
    if (b->core.flag & BAM_FREVERSE) {
        
        alignment.set_sequence(reverse_complement(sequence));
        
        string rev_quality;
        rev_quality.resize(quality.size());
        reverse_copy(quality.begin(), quality.end(), rev_quality.begin());
        alignment.set_quality(rev_quality);
    }
    else {
        
        alignment.set_sequence(sequence);
        alignment.set_quality(quality);
        
    }
    alignment.set_read_paired((b->core.flag & BAM_FPAIRED) != 0);
    
    if (graph != nullptr && bh != nullptr && b->core.tid >= 0) {
        alignment.set_mapping_quality(b->core.qual);
        alignment.set_read_mapped(true);
        // Look for the path handle this is against.
        auto found = tid_path_handle.find(b->core.tid);
        if (found == tid_path_handle.end()) {
            cerr << "[vg::alignment.cpp] error: alignment references path not present in graph: "
                 << bh->target_name[b->core.tid] << endl;
            exit(1);
        }
        mapping_against_path(alignment, b, found->second, graph, b->core.flag & BAM_FREVERSE);
    }
    
    // TODO: htslib doesn't wrap this flag for some reason.
    alignment.set_is_secondary(b->core.flag & BAM_FSECONDARY);
    
    // get the tags
    auto tags = bam_tag_strings(b);
    // handle the tags that are given special fields in GAM
    size_t removed = 0;
    for (size_t i = 0; i < tags.size(); ++i) {
        auto& tag = tags[i];
        auto tag_name = tag.substr(0, 2);
        if (tag_name == "RG") {
            string read_group = tag.substr(5, string::npos);
            alignment.set_read_group(read_group);
            auto it = rg_sample.find(read_group);
            if (it != rg_sample.end()) {
                alignment.set_sample_name(it->second);
            }
            ++removed;
        }
        else if (tag_name == "AS") {
            alignment.set_score(parse<int64_t>(tag.substr(5, string::npos)));
            ++removed;
        }
        else if (removed != 0) {
            tags[i - removed] = std::move(tag);
        }
    }
    
    if (removed != 0) {
        tags.resize(tags.size() - removed);
    }
    
    // save the other tags as an annotation
    if (!tags.empty()) {
        set_annotation(alignment, "tags", tags);
    }

    return alignment;
}

Alignment bam_to_alignment(const bam1_t *b, const map<string, string>& rg_sample, const map<int, path_handle_t>& tid_path_handle) {
    return bam_to_alignment(b, rg_sample, tid_path_handle, nullptr, nullptr);
}

int alignment_to_length(const Alignment& a) {
    int l = 0;
    for (const auto& m : a.path().mapping()) {
        l += to_length(m);
    }
    return l;
}

int alignment_from_length(const Alignment& a) {
    int l = 0;
    for (const auto& m : a.path().mapping()) {
        l += from_length(m);
    }
    return l;
}

Alignment strip_from_start(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    res.set_name(aln.name());
    res.set_score(aln.score());
    res.set_sequence(aln.sequence().substr(drop));
    if (!aln.has_path()) return res;
    *res.mutable_path() = cut_path(aln.path(), drop).second;
    assert(res.has_path());
    if (alignment_to_length(res) != res.sequence().size()) {
        cerr << "failed!!! drop from start 轰" << endl;
        cerr << "drop " << drop << " from start" << endl << pb2json(aln) << endl;
        cerr << "wanted " << aln.sequence().size() - drop << " got " << alignment_to_length(res) << endl;
        cerr << pb2json(res) << endl << endl;
        assert(false);
    }
    return res;
}

Alignment strip_from_end(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    res.set_name(aln.name());
    res.set_score(aln.score());
    //cerr << "drop " << drop << " from end" << endl;
    size_t cut_at = aln.sequence().size()-drop;
    //cerr << "Cut at " << cut_at << endl;
    res.set_sequence(aln.sequence().substr(0, cut_at));
    if (!aln.has_path()) return res;
    *res.mutable_path() = cut_path(aln.path(), cut_at).first;
    assert(res.has_path());
    if (alignment_to_length(res) != res.sequence().size()) {
        cerr << "failed!!! drop from end 轰" << endl;
        cerr << pb2json(res) << endl << endl;
        assert(false);
    }
    return res;
}

Alignment trim_alignment(const Alignment& aln, const Position& pos1, const Position& pos2) {
    // cut the alignment into 3 (possibly empty) pieces
    auto p = cut_path(aln.path(), pos1);
    auto path1 = p.first;
    p = cut_path(p.second, pos2);
    auto path2 = p.first;
    auto path3 = p.second;
    // measure the length of the left and right bits, and use this to trim the current alignment
    auto trimmed = aln;
    if (path1.mapping_size()) {
        trimmed = strip_from_start(trimmed, path_to_length(path1));
    }
    if (path3.mapping_size()) {
        trimmed = strip_from_end(trimmed, path_to_length(path3));
    }
    return trimmed;
}

vector<Alignment> alignment_ends(const Alignment& aln, size_t len1, size_t len2) {
    vector<Alignment> ends;
    ends.push_back(strip_from_end(aln, aln.sequence().size()-len1));
    ends.push_back(strip_from_start(aln, aln.sequence().size()-len2));
    return ends;
}

Alignment alignment_middle(const Alignment& aln, int len) {
    int trim = (aln.sequence().size() - len)/2;
    return strip_from_start(strip_from_end(aln, trim), trim);
}

vector<Alignment> reverse_complement_alignments(const vector<Alignment>& alns, const function<int64_t(int64_t)>& node_length) {
    vector<Alignment> revalns;
    for (auto& aln : alns) {
        revalns.push_back(reverse_complement_alignment(aln, node_length));
    }
    return revalns;
}

Alignment reverse_complement_alignment(const Alignment& aln,
                                       const function<int64_t(id_t)>& node_length) {
    // We're going to reverse the alignment and all its mappings.
    // TODO: should we/can we do this in place?
    
    Alignment reversed = aln;
    reversed.set_sequence(reverse_complement(aln.sequence()));
    string quality = aln.quality();
    std::reverse(quality.begin(), quality.end());
    reversed.set_quality(quality);

    if(aln.has_path()) {
        // Now invert the order of the mappings, and for each mapping, flip the
        // is_reverse flag, and adjust offsets to count from the other end. The
        // edits within mappings also get put in reverse order, and get their
        // sequences reverse complemented.
        *reversed.mutable_path() = reverse_complement_path(aln.path(), node_length);
    }
    
    return reversed;
}
    
void reverse_complement_alignment_in_place(Alignment* aln,
                                           const function<int64_t(id_t)>& node_length) {

    reverse_complement_in_place(*aln->mutable_sequence());
    string* quality = aln->mutable_quality();
    std::reverse(quality->begin(), quality->end());
    
    if (aln->has_path()) {
        reverse_complement_path_in_place(aln->mutable_path(), node_length);
    }
}

// merge that properly handles long indels
// assumes that alignments should line up end-to-end
Alignment merge_alignments(const vector<Alignment>& alns) {

    if (alns.size() == 0) {
        Alignment aln;
        return aln;
    } else if (alns.size() == 1) {
        return alns.front();
    }

    // execute a serial merge
    // buliding up the alignment
    Alignment merged;
    merged.set_name(alns.front().name());

    size_t len = 0;
    for (size_t i = 0; i < alns.size(); ++i) {
        len += alns[i].sequence().size();
    }
    merged.mutable_sequence()->reserve(len);
    if (alns.front().quality().size()) merged.mutable_quality()->reserve(len);

    // get the alignments ready for merge
    for (size_t i = 0; i < alns.size(); ++i) {
        Alignment aln = alns[i];
        if (!aln.has_path()) {
            Mapping m;
            Edit* e = m.add_edit();
            e->set_to_length(aln.sequence().size());
            e->set_sequence(aln.sequence());
            *aln.mutable_path()->add_mapping() = m;
        }
        if (i == 0) {
            merged = aln;
        } else {
            if (!merged.quality().empty()) merged.mutable_quality()->append(aln.quality());
            extend_path(*merged.mutable_path(), aln.path());
            merged.mutable_sequence()->append(aln.sequence());
        }
    }
    return merged;
}

Alignment& extend_alignment(Alignment& a1, const Alignment& a2, bool debug) {
    //if (debug) cerr << "extending alignment " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    a1.set_sequence(a1.sequence() + a2.sequence());
    if (!a1.quality().empty()) a1.set_quality(a1.quality() + a2.quality());
    extend_path(*a1.mutable_path(), a2.path());
    //if (debug) cerr << "extended alignments, result is " << endl << pb2json(a1) << endl;
    return a1;
}

// use a deep copy of the alignments, concatenating them
Alignment merge_alignments(const Alignment& a1, const Alignment& a2, bool debug) {
    //cerr << "overlap is " << overlap << endl;
    // if either doesn't have a path, then treat it like a massive softclip
    if (debug) cerr << "merging alignments " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    // concatenate them
    Alignment a3;
    a3.set_name(a1.name());
    a3.set_sequence(a1.sequence() + a2.sequence());
    *a3.mutable_path() = concat_paths(a1.path(), a2.path());
    if (debug) cerr << "merged alignments, result is " << endl << pb2json(a3) << endl;
    return a3;
}

void translate_nodes(Alignment& a, const unordered_map<id_t, pair<id_t, bool> >& ids, const std::function<size_t(int64_t)>& node_length) {
    Path* path = a.mutable_path();
    for(size_t i = 0; i < path->mapping_size(); i++) {
        // Grab each mapping (includes its position)
        Mapping* mapping = path->mutable_mapping(i);
        auto pos = mapping->position();
        auto oldp = ids.find(pos.node_id());
        if (oldp != ids.end()) {
            auto& old = oldp->second;
            mapping->mutable_position()->set_node_id(old.first);
            if (old.second) {
                mapping->mutable_position()->set_is_reverse(true);
            }
        }
    }
}

void flip_nodes(Alignment& a, const set<int64_t>& ids, const std::function<size_t(int64_t)>& node_length) {
    Path* path = a.mutable_path();
    for(size_t i = 0; i < path->mapping_size(); i++) {
        // Grab each mapping (includes its position)
        Mapping* mapping = path->mutable_mapping(i);
        if(ids.count(mapping->position().node_id())) {
            // We need to flip this mapping
            *mapping = reverse_complement_mapping(*mapping, node_length);
        } 
    }
}

int non_match_start(const Alignment& alignment) {
    int length = 0;
    auto& path = alignment.path();
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                return length;
            }
            length += edit.to_length();
        }
    }
    return length;
}

int non_match_end(const Alignment& alignment) {
    int length = 0;
    auto& path = alignment.path();
    for (int i = path.mapping_size()-1; i >= 0; --i) {
        auto& mapping = path.mapping(i);
        for (int j = mapping.edit_size()-1; j >= 0; --j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                return length;
            }
            length += edit.to_length();
        }
    }
    return length;
}

int softclip_start(const Alignment& alignment) {
    if (alignment.path().mapping_size() > 0) {
        auto& path = alignment.path();
        auto& first_mapping = path.mapping(0);
        auto& first_edit = first_mapping.edit(0);
        if (first_edit.from_length() == 0 && first_edit.to_length() > 0) {
            return first_edit.to_length();
        }
    }
    return 0;
}

int softclip_end(const Alignment& alignment) {
    if (alignment.path().mapping_size() > 0) {
        auto& path = alignment.path();
        auto& last_mapping = path.mapping(path.mapping_size()-1);
        auto& last_edit = last_mapping.edit(last_mapping.edit_size()-1);
        if (last_edit.from_length() == 0 && last_edit.to_length() > 0) {
            return last_edit.to_length();
        }
    }
    return 0;
}

int softclip_trim(Alignment& alignment) {
    // Trim the softclips off of every read
    // Work out were to cut
    int cut_start = softclip_start(alignment);
    int cut_end = softclip_end(alignment);
    // Cut the sequence and quality
    alignment.set_sequence(alignment.sequence().substr(cut_start, alignment.sequence().size() - cut_start - cut_end));
    if (alignment.quality().size() != 0) {
        alignment.set_quality(alignment.quality().substr(cut_start, alignment.quality().size() - cut_start - cut_end));
    }
    // Trim the path
    *alignment.mutable_path() = trim_hanging_ends(alignment.path());
    return cut_start + cut_end;
}

int query_overlap(const Alignment& aln1, const Alignment& aln2) {
    if (!alignment_to_length(aln1) || !alignment_to_length(aln2)
        || !aln1.path().mapping_size() || !aln2.path().mapping_size()
        || aln1.sequence().size() != aln2.sequence().size()) {
        return 0;
    }
    int qb1 = softclip_start(aln1);
    int qe1 = softclip_end(aln1);
    int qb2 = softclip_start(aln2);
    int qe2 = softclip_end(aln2);
    int l = aln1.sequence().size();
    return l - ((qe1 > qe2 ? qe1 : qe2) + (qb1 > qb2 ? qb1 : qb2));
}

int edit_count(const Alignment& alignment) {
    int i = 0;
    auto& path = alignment.path();
    for (int j = path.mapping_size(); j < path.mapping_size(); ++j) {
        i += path.mapping(j).edit_size();
    }
    return i;
}

size_t to_length_after_pos(const Alignment& aln, const Position& pos) {
    return path_to_length(cut_path(aln.path(), pos).second);
}

size_t from_length_after_pos(const Alignment& aln, const Position& pos) {
    return path_from_length(cut_path(aln.path(), pos).second);
}

size_t to_length_before_pos(const Alignment& aln, const Position& pos) {
    return path_to_length(cut_path(aln.path(), pos).first);
}

size_t from_length_before_pos(const Alignment& aln, const Position& pos) {
    return path_from_length(cut_path(aln.path(), pos).first);
}

const string hash_alignment(const Alignment& aln) {
    string data;
    aln.SerializeToString(&data);
    return sha1sum(data);
}

Alignment simplify(const Alignment& a, bool trim_internal_deletions) {
    auto aln = a;
    *aln.mutable_path() = simplify(aln.path(), trim_internal_deletions);
    if (!aln.path().mapping_size()) {
        aln.clear_path();
    }
    return aln;
}
    
void normalize_alignment(Alignment& alignment) {
    
    enum edit_type_t {None, Match, Mismatch, Insert, Delete, N};
    
    size_t cumul_to_length = 0;
    
    // we only build the normalized path if we find things we need to normalize
    // (this makes the whole algorithm a little fucky, but it should be less overhead)
    bool doing_normalization = false;
    Path normalized;
    
    const Path& path = alignment.path();
    const string& seq = alignment.sequence();
    
    auto ensure_init_normalized_path = [&](size_t i, size_t j) {
        // we won't copy the already normalized prefix unless we have to
        if (!doing_normalization) {
            for (size_t k = 0; k < i; k++) {
                *normalized.add_mapping() = path.mapping(k);
            }
            Mapping* mapping = normalized.add_mapping();
            *mapping->mutable_position() = path.mapping(i).position();
            mapping->set_rank(path.mapping_size());
            for (size_t k = 0; k < j; k++) {
                *mapping->add_edit() = path.mapping(i).edit(k);
            }
            doing_normalization = true;
        }
    };
    
    edit_type_t prev = None;
    
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        
        const Mapping& mapping = path.mapping(i);
        prev = None;
        
        if (doing_normalization) {
            // we're maintaining the normalized path, so we need to add mappings
            // as we go
            Mapping* norm_mapping = normalized.add_mapping();
            *norm_mapping->mutable_position() = mapping.position();
            norm_mapping->set_rank(normalized.mapping_size());
        }
        
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            
            const Edit& edit = mapping.edit(j);
            
            if (edit.from_length() > 0 && edit.to_length() == 0) {
                
                if (prev == Delete || doing_normalization) {
                    // we need to modify the normalized path this round
                    ensure_init_normalized_path(i, j);
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (prev == Delete) {
                        // merge with the previous
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                    }
                    else {
                        // just copy
                        *norm_mapping->add_edit() = edit;
                    }
                }
                
                prev = Delete;
            }
            else if (edit.from_length() == 0 && edit.to_length() > 0) {
                
                if (prev == Insert || doing_normalization) {
                    // we need to modify the normalized path this round
                    ensure_init_normalized_path(i, j);
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (prev == Insert) {
                        // merge with the previous
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                        norm_edit->mutable_sequence()->append(edit.sequence());
                    }
                    else {
                        // just copy
                        *norm_mapping->add_edit() = edit;
                    }
                }
                
                cumul_to_length += edit.to_length();
                prev = Insert;
            }
            else {
                auto begin = seq.begin() + cumul_to_length;
                auto end = begin + edit.to_length();
                
                auto first_N = find(begin, end, 'N');
                
                edit_type_t type =  edit.sequence().empty() ? Match : Mismatch;
                
                if (prev == type || first_N != end || doing_normalization) {
                    // we have to do some normalization here
                    ensure_init_normalized_path(i, j);
                    
                    Mapping* norm_mapping = normalized.mutable_mapping(normalized.mapping_size() - 1);
                    if (first_N == end && prev != type) {
                        // just need to copy, no fancy normalization
                        *norm_mapping->add_edit() = edit;
                        prev = type;
                    }
                    else if (first_N == end) {
                        // we need to extend the previous edit, but we don't need
                        // to worry about Ns
                        Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                        norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                        norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                        if (type == Mismatch) {
                            norm_edit->mutable_sequence()->append(edit.sequence());
                        }
                    }
                    else {
                        bool on_Ns = first_N == begin;
                        auto next_pos = begin;
                        // iterate until we've handled the whole edit sequence
                        while (next_pos != end) {
                            // find the next place where we switch from N to non-N or the reverse
                            auto next_end = find_if(next_pos, end, [&](char c) {
                                return c == 'N' != on_Ns;
                            });
                            
                            if ((prev == N && on_Ns) || (prev == type && !on_Ns)) {
                                // we need to merge with the previous edit
                                Edit* norm_edit = norm_mapping->mutable_edit(norm_mapping->edit_size() - 1);
                                norm_edit->set_from_length(norm_edit->from_length() + edit.from_length());
                                norm_edit->set_to_length(norm_edit->to_length() + edit.to_length());
                                
                                // we copy sequence for Ns and for mismatches only
                                if ((prev == N && on_Ns) || (prev == type && !on_Ns && type == Mismatch)) {
                                    norm_edit->mutable_sequence()->append(next_pos, next_end);
                                }
                            }
                            else {
                                // we can just copy
                                Edit* norm_edit = norm_mapping->add_edit();
                                norm_edit->set_from_length(next_end - next_pos);
                                norm_edit->set_to_length(next_end - next_pos);
                                *norm_edit->mutable_sequence() = string(next_pos, next_end);
                            }
                            
                            next_pos = next_end;
                            prev = on_Ns ? N : type;
                            on_Ns = !on_Ns;
                        }
                    }
                }
                else {
                    // no normalization yet
                    prev = type;
                }
                
                cumul_to_length += edit.to_length();
            }
        }
    }
    
    if (doing_normalization) {
        // we found things we needed to normalize away, so we must have built the normalized
        // path, now replace the original with it
        *alignment.mutable_path() = std::move(normalized);
    }
}

bool uses_Us(const Alignment& alignment) {
    
    for (char nt : alignment.sequence()) {
        switch (nt) {
            case 'U':
                return true;
                break;
                
            case 'T':
                return false;
                break;
                
            default:
                break;
        }
    }
    return false;
}

void convert_alignment_char(Alignment& alignment, char from, char to) {
    auto& seq = *alignment.mutable_sequence();
    for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == from) {
            seq[i] = to;
        }
    }
    if (alignment.has_path()) {
        for (Mapping& mapping : *alignment.mutable_path()->mutable_mapping()) {
            for (Edit& edit : *mapping.mutable_edit()) {
                if (!edit.sequence().empty()) {
                    auto& eseq = *edit.mutable_sequence();
                    for (size_t i = 0; i < eseq.size(); ++i) {
                        if (eseq[i] == from) {
                            eseq[i] = to;
                        }
                    }
                }
            }
        }
    }
}

void convert_Us_to_Ts(Alignment& alignment) {
    convert_alignment_char(alignment, 'U', 'T');
}

void convert_Ts_to_Us(Alignment& alignment) {
    convert_alignment_char(alignment, 'T', 'U');
}

map<id_t, int> alignment_quality_per_node(const Alignment& aln) {
    map<id_t, int> quals;
    int to_pos = 0; // offset in quals
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        auto& mapping = aln.path().mapping(i);
        auto to_len = mapping_to_length(mapping);
        if (mapping.has_position()) {
            auto& q = quals[mapping.position().node_id()];
            for (size_t j = 0; j < to_len; ++j) {
                q += aln.quality()[to_pos + j];
            }
        }
        to_pos += mapping_to_length(mapping);
    }
    return quals;
}

string middle_signature(const Alignment& aln, int len) {
    return signature(alignment_middle(aln, len));
}

pair<string, string> middle_signature(const Alignment& aln1, const Alignment& aln2, int len) {
    return make_pair(middle_signature(aln1, len), middle_signature(aln1, len));
}

string signature(const Alignment& aln) {
    stringstream s;
    if (aln.has_path() && aln.path().mapping_size()) {
        auto& pos1 = aln.path().mapping(0).position();
        s << pos1.node_id();
        s << (pos1.is_reverse() ? "-" : "+");
        s << ":" << pos1.offset();
        s << "_";
        auto& last = aln.path().mapping(aln.path().mapping_size()-1);
        auto& pos2 = last.position();
        s << pos2.node_id();
        s << (pos2.is_reverse() ? "-" : "+");
        s << ":" << pos2.offset() + mapping_from_length(last);
    }
    return s.str();
}

pair<string, string> signature(const Alignment& aln1, const Alignment& aln2) {
    return make_pair(signature(aln1), signature(aln2));
}

void parse_bed_regions(istream& bedstream,
                       const PathPositionHandleGraph* graph,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
    parse_bed_regions(bedstream, graph, [&](Alignment& aln) {
        out_alignments->emplace_back(std::move(aln));
    });
}

void parse_bed_regions(istream& bedstream,
                       const PathPositionHandleGraph* graph,
                       const std::function<void(Alignment&)>& callback) {
    
    if (!bedstream) {
        cerr << "Unable to open bed file." << endl;
        return;
    }
    string row;
    string seq;
    // Record start position
    size_t sbuf;
    // Record end position
    size_t ebuf;
    // region information
    string name;
    size_t score = 0;
    string strand;

    // in case we need to look for subpaths, keep the info store for reuse
    unordered_map<string, vector<path_handle_t>> base_path_to_subpaths;
    // to remember which path we've looked for and didn't find (and avoid relooking over and over again)
    unordered_map<string, bool> absent_paths;
    // if the region must be chopped to multiple subranges, use these vector to
    // remember the other path names, start and end positions
    vector<size_t> other_starts;
    vector<size_t> other_ends;
    vector<string> other_seqs;

    for (int line = 1; getline(bedstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
            
        ss >> seq;
        ss >> sbuf;
        ss >> ebuf;

        if (ss.fail()) {
            // Skip lines that can't be parsed
            cerr << "warning: Error parsing bed line " << line << ", skipping: " << row << endl;
            continue;
        } 
        
        if (!graph->has_path(seq)) {
            // This path doesn't exist, and we'll get a segfault or worse if
            // we go look for positions in it.
            // but maybe it's chopped in subranges?
            bool subpath_found = false;
            // first look in our cached subpaths
            // if not there, look in the graph
            if(!base_path_to_subpaths.count(seq) && !absent_paths.count(seq)){
                PathSense sense;
                string sample;
                string locus;
                size_t haplotype;
                size_t phase_block;
                subrange_t subrange;
                PathMetadata::parse_path_name(seq, sense, sample, locus, haplotype, phase_block, subrange);
                if (subrange == PathMetadata::NO_SUBRANGE) {
                    // the path name souldn't describe a subpath
                    graph->for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t& match) {
                        if (graph->get_haplotype(match) != haplotype) {
                            // Skip this haplotype
                            return true;
                        }
                        if (graph->get_phase_block(match) != phase_block) {
                            // Skip this phase block
                            return true;
                        }
                        // we've found a subpath for that base name. save the handle
                        base_path_to_subpaths[seq].push_back(match);
                        return true;
                    });
                }
            }
            if(!base_path_to_subpaths.count(seq)){
                // we've looked for subpaths and couldn't found anything
                cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
                // remember that this path is not in the graph (no need to look it up again)
                absent_paths[seq] = true;
                continue;
            } else {
                // update the path name to the subpaths and adjust sbuf/ebuf based on its offset
                // look for the subpath containing the range [sbuf-ebuf]
                for (auto& path : base_path_to_subpaths[seq]) {
                    subrange_t subrange = graph->get_subrange(path);
                    // the two subrange formats are using different indexing
                    // in path[start-end] start/end are 0-based and the "end" is not included (like the input sbuf/ebuf from BED)
                    // in path[offset] the offset is "1-based"
                    // so to make path[offset] into path[start-end] we need to do path[(offset-1)-(offset-1+path_length)]
                    if (subrange.second == PathMetadata::NO_END_POSITION){
                        if (subrange.first == PathMetadata::NO_END_POSITION){
                            subrange.first = 0;
                        } else {
                            subrange.first--;
                        }
                        subrange.second = subrange.first + graph->get_path_length(path);
                    }
                    // if the subpath overlap with the queried range, save it
                    if(ebuf >= subrange.first && sbuf < subrange.second ){
                        if(sbuf < subrange.first){
                            other_starts.push_back(0);
                        } else {
                            other_starts.push_back(sbuf - subrange.first);
                        }
                        if(ebuf >= subrange.second){
                            other_ends.push_back(subrange.second - subrange.first);
                        } else {
                            other_ends.push_back(ebuf - subrange.first);
                        }
                        other_seqs.push_back(graph->get_path_name(path));
                    }
                }
                if (!other_seqs.empty()){
                    seq = other_seqs.back();
                    other_seqs.pop_back();
                    sbuf = other_starts.back();
                    other_starts.pop_back();
                    ebuf = other_ends.back();
                    other_ends.pop_back();
                } else {
                    // we've looked for overlapping subpaths and couldn't found one
                    cerr << "warning: no overlap found between range " <<
                        sbuf << "-" << ebuf << " and subpaths of \"" <<
                        seq << "\", input line " << line << ", skipping" << endl;
                    continue;
                }
            }
        }

        path_handle_t path_handle = graph->get_path_handle(seq);
        
        if (sbuf >= ebuf && !graph->get_is_circular(path_handle)) {
            // The start of the region can be after the end of the region only if the underlying path is circular.
            // That's not the case, so complain and skip the region.
            cerr << "warning: path \"" << seq << "\" is not circular, skipping end-spanning region on line "
                << line << ": " << row << endl;
            continue;
        }

        if (ebuf > graph->get_path_length(path_handle)) {
            // Could be that the path is chopped but the first subpath is named like the base/full path
            // if so, it was found in the graph but is shorter than the actual full path
            // in case that happened, check for other subpaths that might have that base name
            // if we find one and the range match the queried range, use that
            bool subpath_found = false;
            if(!base_path_to_subpaths.count(seq)){
                // not in our subpath cache, so let's look for subpaths
                PathSense sense;
                string sample;
                string locus;
                size_t haplotype;
                size_t phase_block;
                subrange_t subrange;
                PathMetadata::parse_path_name(seq, sense, sample, locus, haplotype, phase_block, subrange);
                
                if (subrange == PathMetadata::NO_SUBRANGE) {
                    // the path name souldn't describe a subpath
                    graph->for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t& match) {
                        if (graph->get_haplotype(match) != haplotype) {
                            // Skip this haplotype
                            return true;
                        }
                        if (graph->get_phase_block(match) != phase_block) {
                            // Skip this phase block
                            return true;
                        }
                        // we've found a subrange for this base path
                        base_path_to_subpaths[seq].push_back(match);
                        return true;
                    });
                }
            }
            if(!base_path_to_subpaths.count(seq)){
                // Skip ends that are too late
                cerr << "warning: out of range path end " << ebuf << " > " << graph->get_path_length(path_handle)
                     << " in bed line " << line << ", skipping: " << row << endl;
                continue;
            } else {
                // there are subpaths, let's look for the one containing [sbuf-ebuf]
                for (auto& path : base_path_to_subpaths[seq]) {
                    subrange_t subrange = graph->get_subrange(path);
                    // the two subrange formats are using different indexing
                    // in path[start-end] start/end are 0-based and the "end" is not included (like the input sbuf/ebuf from BED)
                    // in path[offset] the offset is "1-based"
                    // so to make path[offset] into path[start-end] we need to do path[(offset-1)-(offset-1+path_length)]
                    if (subrange.second == PathMetadata::NO_END_POSITION){
                        if (subrange.first == PathMetadata::NO_END_POSITION){
                            subrange.first = 0;
                        } else {
                            subrange.first--;
                        }
                        subrange.second = subrange.first + graph->get_path_length(path);
                    }
                    // if the subpath overlap with the queried range, save it
                    if(ebuf >= subrange.first && sbuf < subrange.second ){
                        if(sbuf < subrange.first){
                            other_starts.push_back(0);
                        } else {
                            other_starts.push_back(sbuf - subrange.first);
                        }
                        if(ebuf >= subrange.second){
                            other_ends.push_back(subrange.second - subrange.first);
                        } else {
                            other_ends.push_back(ebuf - subrange.first);
                        }
                        other_seqs.push_back(graph->get_path_name(path));
                    }
                }
                if (!other_seqs.empty()){
                    seq = other_seqs.back();
                    other_seqs.pop_back();
                    path_handle = graph->get_path_handle(seq);
                    sbuf = other_starts.back();
                    other_starts.pop_back();
                    ebuf = other_ends.back();
                    other_ends.pop_back();
                } else {
                    // we've looked for overlapping subpaths and couldn't found one
                    cerr << "warning: no overlap found between range " <<
                        sbuf << "-" << ebuf << " and subpaths of \"" <<
                        seq << "\", input line " << line << ", skipping" << endl;
                    continue;
                }
            }
        }
        
        if (sbuf >= graph->get_path_length(path_handle)) {
            // Skip starts that are too late
            cerr << "warning: out of range path start " << sbuf << " >= " << graph->get_path_length(path_handle)
                << " in bed line " << line << ", skipping: " << row << endl;
            continue;
        }
        
        // Try parsing the optional fields. If they fail, ignore the problem, because they're optional.
        ss >> name;
        ss >> score;
        ss >> strand;

        bool is_reverse = false;
        if(!ss.fail() && strand.compare("-") == 0) {
            is_reverse = true;
        }

        // Make the Alignment
        Alignment alignment = target_alignment(graph, path_handle, sbuf, ebuf, name, is_reverse);
        alignment.set_score(score);
        callback(alignment);

        // if more subpaths need to be written, write them now
        while (!other_seqs.empty()){
            // extract subpath information
            seq = other_seqs.back();
            other_seqs.pop_back();
            sbuf = other_starts.back();
            other_starts.pop_back();
            ebuf = other_ends.back();
            other_ends.pop_back();
            // get path handle and corresponding alignment
            path_handle = graph->get_path_handle(seq);
            alignment = target_alignment(graph, path_handle, sbuf, ebuf, name, is_reverse);
            alignment.set_score(score);
            callback(alignment);
        }
    }
}

void parse_gff_regions(istream& gffstream,
                       const PathPositionHandleGraph* graph,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
    parse_gff_regions(gffstream, graph, [&](Alignment& aln) {
        out_alignments->emplace_back(std::move(aln));
    });
}

void parse_gff_regions(istream& gffstream,
                       const PathPositionHandleGraph* graph,
                       const std::function<void(Alignment&)>& callback) {


    if (!gffstream) {
        cerr << "Unable to open gff3/gtf file." << endl;
        return;
    }
    string row;
    string seq;
    string source;
    string type;
    string buf;
    size_t sbuf;
    size_t ebuf;
    string name = "";
    string score;
    string strand;
    string num;
    string annotations;

    // in case we need to look for subpaths, keep the info store for reuse
    unordered_map<string, vector<path_handle_t>> base_path_to_subpaths;
    // to remember which path we've looked for and didn't find (and avoid relooking over and over again)
    unordered_map<string, bool> absent_paths;
    // if the region must be chopped to multiple subranges, use these vector to
    // remember the other path names, start and end positions
    vector<size_t> other_starts;
    vector<size_t> other_ends;
    vector<string> other_seqs;

    for (int line = 1; getline(gffstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        getline(ss, seq, '\t');
        getline(ss, source, '\t');
        getline(ss, type, '\t');
        getline(ss, buf, '\t');
        // Convert to 0-based 
        sbuf = atoi(buf.c_str()) - 1;
        getline(ss, buf, '\t');
        // 1-based inclusive == 0-based exclusive
        ebuf = atoi(buf.c_str());

        if (ss.fail() || !(sbuf < ebuf)) {
            cerr << "Error parsing gtf/gff line " << line << ": " << row << endl;
            continue;
        }
        
        getline(ss, score, '\t');
        getline(ss, strand, '\t');
        getline(ss, num, '\t');
        getline(ss, annotations, '\t');

        // look for the "Name" info
        vector<string> vals = split(annotations, ";");
        string name = "";
        for (auto& s : vals) {
            if (s.find("Name=") == 0) {
                name = s.substr(5);
            }
        }

        // Skips annotations where the name can not be parsed. Empty names can 
        // results in undefinable behavior downstream. 
        if (name.empty()) {
            cerr << "warning: could not parse annotation name (Name=), skipping line " << line << endl;  
            continue;              
        }

        bool is_reverse = false;
        if(!ss.fail() && strand.compare("-") == 0) {
            is_reverse = true;
        }

        if (!graph->has_path(seq)) {
            // This path doesn't exist, and we'll get a segfault or worse if
            // we go look for positions in it.
            // but maybe it's chopped in subranges?
            bool subpath_found = false;
            // first look in our cached subpaths
            // if not there, look in the graph
            if(!base_path_to_subpaths.count(seq) && !absent_paths.count(seq)){
                PathSense sense;
                string sample;
                string locus;
                size_t haplotype;
                size_t phase_block;
                subrange_t subrange;
                PathMetadata::parse_path_name(seq, sense, sample, locus, haplotype, phase_block, subrange);
                if (subrange == PathMetadata::NO_SUBRANGE) {
                    // the path name souldn't describe a subpath
                    graph->for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t& match) {
                        if (graph->get_haplotype(match) != haplotype) {
                            // Skip this haplotype
                            return true;
                        }
                        if (graph->get_phase_block(match) != phase_block) {
                            // Skip this phase block
                            return true;
                        }
                        // we've found a subpath for that base name. save the handle
                        base_path_to_subpaths[seq].push_back(match);
                        return true;
                    });
                }
            }
            if(!base_path_to_subpaths.count(seq)){
                // we've looked for subpaths and couldn't found anything
                cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
                // remember that this path is not in the graph (no need to look it up again)
                absent_paths[seq] = true;
                continue;
            } else {
                // update the path name to the subpaths and adjust sbuf/ebuf based on its offset
                // look for the subpath containing the range [sbuf-ebuf]
                for (auto& path : base_path_to_subpaths[seq]) {
                    subrange_t subrange = graph->get_subrange(path);
                    // the two subrange formats are using different indexing
                    // in path[start-end] start/end are 0-based and the "end" is not included (like the input sbuf/ebuf from BED)
                    // in path[offset] the offset is "1-based"
                    // so to make path[offset] into path[start-end] we need to do path[(offset-1)-(offset-1+path_length)]
                    if (subrange.second == PathMetadata::NO_END_POSITION){
                        if (subrange.first == PathMetadata::NO_END_POSITION){
                            subrange.first = 0;
                        } else {
                            subrange.first--;
                        }
                        subrange.second = subrange.first + graph->get_path_length(path);
                    }
                    // if the subpath overlap with the queried range, save it
                    if(ebuf >= subrange.first && sbuf < subrange.second ){
                        if(sbuf < subrange.first){
                            other_starts.push_back(0);
                        } else {
                            other_starts.push_back(sbuf - subrange.first);
                        }
                        if(ebuf >= subrange.second){
                            other_ends.push_back(subrange.second - subrange.first);
                        } else {
                            other_ends.push_back(ebuf - subrange.first);
                        }
                        other_seqs.push_back(graph->get_path_name(path));
                    }
                }
                if (!other_seqs.empty()){
                    seq = other_seqs.back();
                    other_seqs.pop_back();
                    sbuf = other_starts.back();
                    other_starts.pop_back();
                    ebuf = other_ends.back();
                    other_ends.pop_back();
                } else {
                    // we've looked for overlapping subpaths and couldn't found one
                    cerr << "warning: no overlap found between range " <<
                        sbuf << "-" << ebuf << " and subpaths of \"" <<
                        seq << "\", input line " << line << ", skipping" << endl;
                    continue;
                }
            }
        }

        path_handle_t path_handle = graph->get_path_handle(seq);
        
        if (ebuf > graph->get_path_length(path_handle)) {
            // Could be that the path is chopped but the first subpath is named like the base/full path
            // if so, it was found in the graph but is shorter than the actual full path
            // in case that happened, check for other subpaths that might have that base name
            // if we find one and the range match the queried range, use that
            bool subpath_found = false;
            if(!base_path_to_subpaths.count(seq)){
                // not in our subpath cache, so let's look for subpaths
                PathSense sense;
                string sample;
                string locus;
                size_t haplotype;
                size_t phase_block;
                subrange_t subrange;
                PathMetadata::parse_path_name(seq, sense, sample, locus, haplotype, phase_block, subrange);
                
                if (subrange == PathMetadata::NO_SUBRANGE) {
                    // the path name souldn't describe a subpath
                    graph->for_each_path_matching({sense}, {sample}, {locus}, [&](const path_handle_t& match) {
                        if (graph->get_haplotype(match) != haplotype) {
                            // Skip this haplotype
                            return true;
                        }
                        if (graph->get_phase_block(match) != phase_block) {
                            // Skip this phase block
                            return true;
                        }
                        // we've found a subrange for this base path
                        base_path_to_subpaths[seq].push_back(match);
                        return true;
                    });
                }
            }
            if(!base_path_to_subpaths.count(seq)){
                // Skip ends that are too late
                cerr << "warning: out of range path end " << ebuf << " > " << graph->get_path_length(path_handle)
                     << " in bed line " << line << ", skipping: " << row << endl;
                continue;
            } else {
                // there are subpaths, let's look for the one containing [sbuf-ebuf]
                for (auto& path : base_path_to_subpaths[seq]) {
                    subrange_t subrange = graph->get_subrange(path);
                    // the two subrange formats are using different indexing
                    // in path[start-end] start/end are 0-based and the "end" is not included (like the input sbuf/ebuf from BED)
                    // in path[offset] the offset is "1-based"
                    // so to make path[offset] into path[start-end] we need to do path[(offset-1)-(offset-1+path_length)]
                    if (subrange.second == PathMetadata::NO_END_POSITION){
                        if (subrange.first == PathMetadata::NO_END_POSITION){
                            subrange.first = 0;
                        } else {
                            subrange.first--;
                        }
                        subrange.second = subrange.first + graph->get_path_length(path);
                    }
                    // if the subpath overlap with the queried range, save it
                    if(ebuf >= subrange.first && sbuf < subrange.second ){
                        if(sbuf < subrange.first){
                            other_starts.push_back(0);
                        } else {
                            other_starts.push_back(sbuf - subrange.first);
                        }
                        if(ebuf >= subrange.second){
                            other_ends.push_back(subrange.second - subrange.first);
                        } else {
                            other_ends.push_back(ebuf - subrange.first);
                        }
                        other_seqs.push_back(graph->get_path_name(path));
                    }
                }
                if (!other_seqs.empty()){
                    seq = other_seqs.back();
                    other_seqs.pop_back();
                    path_handle = graph->get_path_handle(seq);
                    sbuf = other_starts.back();
                    other_starts.pop_back();
                    ebuf = other_ends.back();
                    other_ends.pop_back();
                } else {
                    // we've looked for overlapping subpaths and couldn't found one
                    cerr << "warning: no overlap found between range " <<
                        sbuf << "-" << ebuf << " and subpaths of \"" <<
                        seq << "\", input line " << line << ", skipping" << endl;
                    continue;
                }
            }
        }
        
        if (sbuf >= graph->get_path_length(path_handle)) {
            // Skip starts that are too late
            cerr << "warning: out of range path start " << sbuf << " >= " << graph->get_path_length(path_handle)
                << " in bed line " << line << ", skipping: " << row << endl;
            continue;
        }
        
        Alignment alignment = target_alignment(graph, graph->get_path_handle(seq), sbuf, ebuf, name, is_reverse);
        callback(alignment);

        // if more subpaths need to be written, write them now
        while (!other_seqs.empty()){
            // extract subpath information
            seq = other_seqs.back();
            other_seqs.pop_back();
            sbuf = other_starts.back();
            other_starts.pop_back();
            ebuf = other_ends.back();
            other_ends.pop_back();
            // get alignment
            alignment = target_alignment(graph, graph->get_path_handle(seq), sbuf, ebuf, name, is_reverse);
            callback(alignment);
        }
    }
}

Position alignment_start(const Alignment& aln) {
    Position pos;
    if (aln.path().mapping_size()) {
        pos = aln.path().mapping(0).position();
    }
    return pos;
}

Position alignment_end(const Alignment& aln) {
    Position pos;
    if (aln.path().mapping_size()) {
        auto& last = aln.path().mapping(aln.path().mapping_size()-1);
        pos = last.position();
        pos.set_offset(pos.offset() + mapping_from_length(last));
    }
    return pos;
}

map<string ,vector<pair<size_t, bool> > > alignment_refpos_to_path_offsets(const Alignment& aln) {
    map<string, vector<pair<size_t, bool> > > offsets;
    for (auto& refpos : aln.refpos()) {
        offsets[refpos.name()].push_back(make_pair(refpos.offset(), refpos.is_reverse()));
    }
    return offsets;
}

void alignment_set_distance_to_correct(Alignment& aln, const Alignment& base, const unordered_map<string, string>* translation) {
    auto base_offsets = alignment_refpos_to_path_offsets(base);
    return alignment_set_distance_to_correct(aln, base_offsets, translation);
}

void alignment_set_distance_to_correct(Alignment& aln, const map<string ,vector<pair<size_t, bool> > >& base_offsets, const unordered_map<string, string>* translation) {
    auto aln_offsets = alignment_refpos_to_path_offsets(aln);
    // bail out if we can't compare
    if (!(aln_offsets.size() && base_offsets.size())) return;
    // otherwise find the minimum distance and relative orientation
    Position result;
    size_t min_distance = std::numeric_limits<size_t>::max();
    for (auto& path : aln_offsets) {
        auto name = path.first;
        if (translation) {
            // See if we need to translate the name of the path
            auto found = translation->find(name);
            if (found != translation->end()) {
                // We have a replacement so apply it.
                name = found->second;
            }
        }
        auto& aln_positions = path.second;
        auto f = base_offsets.find(name);
        if (f == base_offsets.end()) continue;
        auto& base_positions = f->second;
        for (auto& p1 : aln_positions) {
            for (auto& p2 : base_positions) {
                // disable relative inversions
                if (p1.second != p2.second) continue;
                // are they in the same orientation?
                size_t dist = abs((int64_t)p1.first - (int64_t)p2.first);
                if (dist < min_distance) {
                    min_distance = dist;
                    result.set_name(name);
                    result.set_is_reverse(p1.second != p2.second);
                    result.set_offset(dist);
                }
            }
        }
    }
    // set the distance to correct if we got one
    if (min_distance < std::numeric_limits<size_t>::max()) {
        *aln.mutable_to_correct() = result;
    }
}

AlignmentValidity alignment_is_valid(const Alignment& aln, const HandleGraph* hgraph, bool check_sequence) {
    size_t read_idx = 0;
    for (size_t i = 0; i < aln.path().mapping_size(); ++i) {
        const Mapping& mapping = aln.path().mapping(i);
        if (!hgraph->has_node(mapping.position().node_id())) {
            std::stringstream ss;
            ss << "Node " << mapping.position().node_id() << " not found in graph";
            return {
                AlignmentValidity::NODE_MISSING,
                i,
                ss.str()
            };
        }
        size_t node_len = hgraph->get_length(hgraph->get_handle(mapping.position().node_id()));
        if (mapping_from_length(mapping) + mapping.position().offset() > node_len) {
            std::stringstream ss;
            ss << "Length of node "
               << mapping.position().node_id() << " (" << node_len << ") exceeded by Mapping with offset "
               << mapping.position().offset() << " and from-length " << mapping_from_length(mapping);
            return {
                AlignmentValidity::NODE_TOO_SHORT,
                i,
                ss.str()
            };
        }
        if (check_sequence) {
            size_t node_idx = mapping.position().offset();
            auto node_seq = hgraph->get_sequence(hgraph->get_handle(mapping.position().node_id(),
                                                                    mapping.position().is_reverse()));
            for (size_t j = 0; j < mapping.edit_size(); ++j) {
                const auto& edit = mapping.edit(j);
                if (edit.to_length() == edit.from_length() && edit.from_length() != 0) {
                    assert(edit.sequence().size() == edit.to_length() || edit.sequence().empty());
                    for (size_t k = 0; k < edit.to_length(); ++k) {
                        // check match/mismatch state between read and ref
                        if ((aln.sequence()[read_idx + k] == node_seq[node_idx + k]) != edit.sequence().empty()) {
                            std::stringstream ss;
                            ss << "Edit erroneously claims " << (edit.sequence().empty() ? "match" : "mismatch") << " on node " << mapping.position().node_id() << " between node position " << (node_idx + k) << " and edit " << j << ", position " << k << " on " << (mapping.position().is_reverse() ? "reverse" : "forward") << " strand";
                            return {
                                AlignmentValidity::SEQ_DOES_NOT_MATCH,
                                i,
                                ss.str()
                            };
                        }
                        if (!edit.sequence().empty() && edit.sequence()[k] != aln.sequence()[read_idx + k]) {
                            // compare mismatched sequence to the read
                            std::stringstream ss;
                            ss << "Edit sequence (" << edit.sequence() << ") at position " << k << " does not match read sequence (" << aln.sequence() << ") at position " << (read_idx + k);
                            return {
                                AlignmentValidity::SEQ_DOES_NOT_MATCH,
                                i,
                                ss.str()
                            };
                        }
                    }
                }
                else if (edit.from_length() == 0 && edit.to_length() != 0) {
                    // compare inserted sequence to read
                    assert(edit.sequence().size() == edit.to_length());
                    for (size_t k = 0; k < edit.to_length(); ++k) {
                        if (edit.sequence()[k] != aln.sequence()[read_idx + k]) {
                            std::stringstream ss;
                            ss << "Read sequence (" << aln.sequence() << ") at position " << (read_idx + k) << " does not match insert sequence of edit (" << edit.sequence() << ") at position " << k;
                            return {
                                AlignmentValidity::SEQ_DOES_NOT_MATCH,
                                i,
                                ss.str()
                            };
                        }
                    }
                }
                else {
                    assert(edit.from_length() != 0 && edit.to_length() == 0);
                }
                
                node_idx += edit.from_length();
                read_idx += edit.to_length();
            }
        }
    }
    return {AlignmentValidity::OK};
}

Alignment target_alignment(const PathPositionHandleGraph* graph, const path_handle_t& path, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse, Mapping& cigar_mapping) {
    Alignment aln;
    
    // How long is the path?
    auto path_len = graph->get_path_length(path);
    
    if (pos2 < pos1) {
        // Looks like we want to span the origin of a circular path
        if (!graph->get_is_circular(path)) {
            // But the path isn't circular, which is a problem
            throw runtime_error("Cannot extract Alignment from " + to_string(pos1) +
                                " to " + to_string(pos2) + " across the junction of non-circular path " +
                                graph->get_path_name(path));
        }
        
        if (pos1 >= path_len) {
            // We want to start off the end of the path, which is no good.
            throw runtime_error("Cannot extract Alignment starting at " + to_string(pos1) +
                                " which is past end " + to_string(path_len) + " of path " +
                                graph->get_path_name(path));
        }
        
        if (pos2 > path_len) {
            // We want to end off the end of the path, which is no good either.
            throw runtime_error("Cannot extract Alignment ending at " + to_string(pos2) +
                                " which is past end " + to_string(path_len) + " of path " +
                                graph->get_path_name(path));
        }
        
        // Split the provided Mapping of edits at the path end/start junction
        auto part_mappings = cut_mapping_offset(cigar_mapping, path_len - pos1);
        
        // We extract from pos1 to the end
        Alignment aln1 = target_alignment(graph, path, pos1, path_len, feature, is_reverse, part_mappings.first);
        
        // And then from the start to pos2
        Alignment aln2 = target_alignment(graph, path, 0, pos2, feature, is_reverse, part_mappings.second);
        
        if (is_reverse) {
            // The alignments were flipped, so the second has to be first
            return merge_alignments(aln2, aln1);
        } else {
            // The alignments get merged in the same order
            return merge_alignments(aln1, aln2);
        }
    }
    
    // Otherwise, the base case is that we don't go over the circular path junction
    
    if (pos1 >= path_len) {
        throw runtime_error("Cannot extract Alignment starting at " + to_string(pos1) +
                            " which is past end " + to_string(path_len) + " of path " +
                            graph->get_path_name(path));
    }
    if (pos2 > path_len) {
        throw runtime_error("Cannot extract Alignment ending at " + to_string(pos2) +
                            " which is past end " + to_string(path_len) + " of path " +
                            graph->get_path_name(path));
    }
    
    step_handle_t step = graph->get_step_at_position(path, pos1);
    
    size_t edit_idx = 0;
    size_t offset_in_edit = 0;
    size_t node_pos = pos1 - graph->get_position_of_step(step);
    while (edit_idx < cigar_mapping.edit_size()) {
        if (step == graph->path_end(path)) {
            const auto& edit = cigar_mapping.edit(edit_idx);
            if (edit.from_length() == 0 && aln.path().mapping_size() != 0) {
                // This is a softclip off the end of the contig.
                // We can add it to the last mapping as an edit
                assert(offset_in_edit == 0);
                Mapping* last_mapping = aln.mutable_path()->mutable_mapping(aln.path().mapping_size() - 1);
                *last_mapping->add_edit() = edit;
                ++edit_idx;
                continue;
            } else {
                // We've gone off the end of the contig with something other than a softclip
                throw std::runtime_error("Reached unexpected end of path " + graph->get_path_name(path) +
                                         " at edit " + std::to_string(edit_idx) +
                                         "/" + std::to_string(cigar_mapping.edit_size()) +
                                         " for alignment of feature " + feature);
            }
        }
        handle_t h = graph->get_handle_of_step(step);
        string seq = graph->get_sequence(h);

        auto mapping = aln.mutable_path()->add_mapping();

        mapping->mutable_position()->set_node_id(graph->get_id(h));
        mapping->mutable_position()->set_is_reverse(graph->get_is_reverse(h));
        mapping->mutable_position()->set_offset(node_pos);
        mapping->set_rank(aln.path().mapping_size());

        while (edit_idx < cigar_mapping.edit_size() && node_pos < seq.size()) {

            const auto& edit = cigar_mapping.edit(edit_idx);

            if (edit.from_length() == edit.to_length()) {
                // match/mismatch -- need to check

                // end at the sooner of 1) the end of the edit and 2) the end of the node
                size_t node_aln_len = min<size_t>(edit.from_length() - offset_in_edit, seq.size() - node_pos);

                // iterate through node and edit up to the limit
                Edit* new_edit = nullptr;
                for (size_t i = 0; i < node_aln_len; ++i, ++offset_in_edit, ++node_pos) {

                    bool match = (edit.sequence()[offset_in_edit] == seq[node_pos]);
                    if (!new_edit || match != new_edit->sequence().empty()) {
                        // current edit is of the wrong type or doesn't exist
                        new_edit = mapping->add_edit();
                    }
                    new_edit->set_from_length(new_edit->from_length() + 1);
                    new_edit->set_to_length(new_edit->to_length() + 1);
                    if (!match) {
                        new_edit->mutable_sequence()->push_back(edit.sequence()[offset_in_edit]);
                    }
                }
                if (offset_in_edit == edit.from_length()) {
                    ++edit_idx;
                    offset_in_edit = 0;
                }
            }
            else if (edit.from_length() == 0) {
                // insertion
                assert(offset_in_edit == 0);
                *mapping->add_edit() = edit;
                ++edit_idx;
            }
            else {
                // deletion
                auto new_edit = mapping->add_edit();
                size_t edit_remaining = edit.from_length() - offset_in_edit;
                size_t node_remaining = seq.size() - node_pos;
                if (edit_remaining <= node_remaining) {
                    // we hit the end of the edit before the end of the node
                    new_edit->set_from_length(edit_remaining);
                    ++edit_idx;
                    offset_in_edit = 0;
                }
                else {
                    // we hit the end of the node before the end of the edit
                    new_edit->set_from_length(node_remaining);
                    offset_in_edit += new_edit->from_length();
                }
                node_pos += new_edit->from_length();
            }
        }

        step = graph->get_next_step(step);
        node_pos = 0;
    }
    
    aln.set_name(feature);
    if (is_reverse) {
        reverse_complement_alignment_in_place(&aln, [&](vg::id_t node_id) { return graph->get_length(graph->get_handle(node_id)); });
    }
    return aln;
}

Alignment target_alignment(const PathPositionHandleGraph* graph, const path_handle_t& path, size_t pos1, size_t pos2,
                           const string& feature, bool is_reverse) {
    Alignment aln;

    
    if (pos2 < pos1) {
        // Looks like we want to span the origin of a circular path
        if (!graph->get_is_circular(path)) {
            // But the path isn't circular, which is a problem
            throw runtime_error("Cannot extract Alignment from " + to_string(pos1) +
                                " to " + to_string(pos2) + " across the junction of non-circular path " +
                                graph->get_path_name(path));
        }
        
        // How long is the path?
        auto path_len = graph->get_path_length(path);
        
        if (pos1 >= path_len) {
            // We want to start off the end of the path, which is no good.
            throw runtime_error("Cannot extract Alignment starting at " + to_string(pos1) +
                                " which is past end " + to_string(path_len) + " of path " +
                                graph->get_path_name(path));
        }
        
        if (pos2 > path_len) {
            // We want to end off the end of the path, which is no good either.
            throw runtime_error("Cannot extract Alignment ending at " + to_string(pos2) +
                                " which is past end " + to_string(path_len) + " of path " +
                                graph->get_path_name(path));
        }
        
        // We extract from pos1 to the end
        Alignment aln1 = target_alignment(graph, path, pos1, path_len, feature, is_reverse);
        
        // And then from the start to pos2
        Alignment aln2 = target_alignment(graph, path, 0, pos2, feature, is_reverse);
        
        if (is_reverse) {
            // The alignments were flipped, so the second has to be first
            return merge_alignments(aln2, aln1);
        } else {
            // The alignments get merged in the same order
            return merge_alignments(aln1, aln2);
        }
    }
    
    // If we get here, we do the normal non-circular path case.
    
    step_handle_t step = graph->get_step_at_position(path, pos1);
    size_t step_start = graph->get_position_of_step(step);
    handle_t handle = graph->get_handle_of_step(step);
    
    int64_t trim_start = pos1 - step_start;
    {
        Mapping* first_mapping = aln.mutable_path()->add_mapping();
        first_mapping->mutable_position()->set_node_id(graph->get_id(handle));
        first_mapping->mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        first_mapping->mutable_position()->set_offset(trim_start);
        
        Edit* e = first_mapping->add_edit();
        size_t edit_len = min<size_t>(graph->get_length(handle) - trim_start, pos2 - pos1);
        e->set_from_length(edit_len);
        e->set_to_length(edit_len);
    }
    // get p to point to the next step (or past it, if we're a feature on a single node)
    int64_t p = step_start + graph->get_length(handle);
    step = graph->get_next_step(step);
    while (p < pos2) {
        handle = graph->get_handle_of_step(step);
        
        Mapping* m = aln.mutable_path()->add_mapping();
        m->mutable_position()->set_node_id(graph->get_id(handle));
        m->mutable_position()->set_is_reverse(graph->get_is_reverse(handle));
        
        Edit* e = m->add_edit();
        size_t edit_len = min<size_t>(graph->get_length(handle), pos2 - p);
        e->set_from_length(edit_len);
        e->set_to_length(edit_len);
        
        p += graph->get_length(handle);
        step = graph->get_next_step(step);
    }
    
    aln.set_name(feature);
    if (is_reverse) {
        reverse_complement_alignment_in_place(&aln, [&](vg::id_t node_id) { return graph->get_length(graph->get_handle(node_id)); });
    }
    return aln;
}
}
