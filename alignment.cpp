#include "alignment.hpp"
#include "stream.hpp"

namespace vg {

int hts_for_each(string& filename, function<void(Alignment&)> lambda) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        Alignment a = bam_to_alignment(b, rg_sample);
        lambda(a);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

}

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);

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
#pragma omp critical (hts_input)
            if (more_data) {
                more_data = sam_read1(in, hdr, b) >= 0;
            }
            if (more_data) {
                Alignment a = bam_to_alignment(b, rg_sample);
                lambda(a);
            }
        }
    }

    for (auto& b : bs) bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    return 1;

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
                             map<string, int64_t>& path_length,
                             map<string, string>& rg_sample) {
    stringstream hdr;
    hdr << "@HD\tVN:1.5\tSO:unknown\n";
    for (auto& p : path_length) {
        hdr << "@SQ\tSN:" << p.first << "\t" << "LN:" << p.second << "\n";
    }
    for (auto& s : rg_sample) {
        hdr << "@RG\tID:" << s.first << "\t" << "SM:" << s.second << "\n";
    }
    hdr << "@PG\tID:0\tPN:vg\n";
    header = hdr.str();
    string sam = "data:" + header;
    samFile *in = sam_open(sam.c_str(), "r");
    bam_hdr_t *h = sam_hdr_read(in);
    sam_close(in);
    return h;
}

bool get_next_alignment_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& alignment) {

    alignment.Clear();

    // handle name
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        string name = buffer;
        name = name.substr(1); // trim off leading @
        // XXX todo trim trailing /1 /2
        alignment.set_name(name);
    } else { return false; }
    // handle sequence
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        alignment.set_sequence(buffer);
    } else {
        cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
    }
    // handle "+" sep
    if (0!=gzgets(fp,buffer,len)) {
    } else {
        cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
    }
    // handle quality
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        string quality = string_quality_char_to_short(buffer);
        alignment.set_quality(quality);
    } else {
        cerr << "[vg::alignment.cpp] error: incomplete fastq record" << endl; exit(1);
    }

    return true;

}

bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp, buffer, len, mate1) && get_next_alignment_from_fastq(fp, buffer, len, mate2);
}

bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp1, buffer, len, mate1) && get_next_alignment_from_fastq(fp2, buffer, len, mate2);
}


size_t fastq_unpaired_for_each_parallel(string& filename, function<void(Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    int thread_count = get_thread_count();
    //vector<Alignment> alns; alns.resize(thread_count);
    vector<char*> bufs; bufs.resize(thread_count);
    for (auto& buf : bufs) {
        buf = new char[len];
    }
    bool more_data = true;
#pragma omp parallel shared(fp, more_data, bufs)
    {
        int tid = omp_get_thread_num();
        while (more_data) {
            Alignment aln;
            char* buf = bufs[tid];
#pragma omp critical (fastq_input)
            if (more_data) {
                more_data = get_next_alignment_from_fastq(fp, buf, len, aln);
                nLines++;
            }
            if (more_data) {
                lambda(aln);
            }
        }
    }
    for (auto& buf : bufs) {
        delete buf;
    }
    gzclose(fp);
    return nLines;
}

size_t fastq_paired_interleaved_for_each_parallel(string& filename, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    int thread_count = get_thread_count();
    //vector<Alignment> alns; alns.resize(thread_count);
    vector<char*> bufs; bufs.resize(thread_count);
    for (auto& buf : bufs) {
        buf = new char[len];
    }
    bool more_data = true;
#pragma omp parallel shared(fp, more_data, bufs)
    {
        int tid = omp_get_thread_num();
        while (more_data) {
            Alignment mate1, mate2;
            char* buf = bufs[tid];
#pragma omp critical (fastq_input)
            if (more_data) {
                more_data = get_next_interleaved_alignment_pair_from_fastq(fp, buf, len, mate1, mate2);
                nLines++;
            }
            if (more_data) {
                lambda(mate1, mate2);
            }
        }
    }
    for (auto& buf : bufs) {
        delete buf;
    }
    gzclose(fp);
    return nLines;
}

size_t fastq_paired_two_files_for_each_parallel(string& file1, string& file2, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    int thread_count = get_thread_count();
    //vector<Alignment> alns; alns.resize(thread_count);
    vector<char*> bufs; bufs.resize(thread_count);
    for (auto& buf : bufs) {
        buf = new char[len];
    }
    bool more_data = true;
#pragma omp parallel shared(fp1, fp2, more_data, bufs)
    {
        int tid = omp_get_thread_num();
        while (more_data) {
            Alignment mate1, mate2;
            char* buf = bufs[tid];
#pragma omp critical (fastq_input)
            if (more_data) {
                more_data = get_next_alignment_pair_from_fastqs(fp1, fp2, buf, len, mate1, mate2);
                nLines++;
            }
            if (more_data) {
                lambda(mate1, mate2);
            }
        }
    }
    for (auto& buf : bufs) {
        delete buf;
    }
    gzclose(fp1);
    gzclose(fp2);
    return nLines;
}



size_t fastq_unpaired_for_each(string& filename, function<void(Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment alignment;
    while(get_next_alignment_from_fastq(fp, buffer, len, alignment)) {
        lambda(alignment);
        nLines++;
    }
    gzclose(fp);
    delete buffer;
    return nLines;
}

size_t fastq_paired_interleaved_for_each(string& filename, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    size_t len = 2 << 18; // 256k
    size_t nLines = 0;
    char *buffer = new char[len];
    Alignment mate1, mate2;
    while(get_next_interleaved_alignment_pair_from_fastq(fp, buffer, len, mate1, mate2)) {
        lambda(mate1, mate2);
        nLines++;
    }
    gzclose(fp);
    delete buffer;
    return nLines;
}

size_t fastq_paired_two_files_for_each(string& file1, string& file2, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp1 = (file1 != "-") ? gzopen(file1.c_str(), "r") : gzdopen(fileno(stdin), "r");
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
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
    delete buffer;
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

void write_alignments(std::ostream& out, vector<Alignment>& buf) {
    function<Alignment(uint64_t)> lambda =
        [&buf] (uint64_t n) {
        return buf[n];
    };
    stream::write(cout, buf.size(), lambda);
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

// remember to clean up with bam_destroy1(b);
bam1_t* alignment_to_bam(const string& sam_header,
                         const Alignment& alignment,
                         const string& refseq,
                         const int32_t refpos,
                         const string& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         const int32_t tlen) {

    assert(!sam_header.empty());
    string sam_file = "data:" + sam_header + alignment_to_sam(alignment, refseq, refpos, cigar, mateseq, matepos, tlen);
    const char* sam = sam_file.c_str();
    samFile *in = sam_open(sam, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    if (sam_read1(in, header, aln) >= 0) {
        bam_hdr_destroy(header);
        sam_close(in); // clean up
        return aln;
    } else {
        cerr << "[vg::alignment] Failure to parse SAM record" << endl
             << sam << endl;
        exit(1);
    }
}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen) {
    stringstream sam;

    sam << (!alignment.name().empty() ? alignment.name() : "*") << "\t"
        << sam_flag(alignment) << "\t"
        << (refseq.empty() ? "*" : refseq) << "\t"
        << refpos + 1 << "\t"
        //<< (alignment.path().mapping_size() ? refpos + 1 : 0) << "\t" // positions are 1-based in SAM, 0 means unmapped
        << alignment.mapping_quality() << "\t"
        << (alignment.has_path() && alignment.path().mapping_size() ? cigar : "*") << "\t"
        << (mateseq == refseq ? "=" : mateseq) << "\t"
        << matepos + 1 << "\t"
        << tlen << "\t"
        << (!alignment.sequence().empty() ? alignment.sequence() : "*") << "\t";
    // hack much?
    if (!alignment.quality().empty()) {
        const string& quality = alignment.quality();
        for (int i = 0; i < quality.size(); ++i) {
            sam << quality_short_to_char(quality[i]);
        }
    } else {
        sam << "*";
        //sam << string(alignment.sequence().size(), 'I');
    }
    //<< (alignment.has_quality() ? string_quality_short_to_char(alignment.quality()) : string(alignment.sequence().size(), 'I'));
    if (!alignment.read_group().empty()) sam << "\tRG:Z:" << alignment.read_group();
    sam << "\n";
    return sam.str();
}

// act like the path this is against is the reference
// and generate an equivalent cigar
string cigar_against_path(const Alignment& alignment) {
    vector<pair<int, char> > cigar;
    if (!alignment.has_path()) return "";
    const Path& path = alignment.path();
    int l = 0;
    for (const auto& mapping : path.mapping()) {
        mapping_cigar(mapping, cigar);
    }
    return cigar_string(cigar);
}

int32_t sam_flag(const Alignment& alignment) {
    int16_t flag = 0;

    if (alignment.score() == 0) {
        // unmapped
        flag |= BAM_FUNMAP;
    } else {
        // correctly aligned
        flag |= BAM_FPROPER_PAIR;
    }
    if (alignment.is_reverse()) {
        flag |= BAM_FREVERSE;
    }
    return flag;
}

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample) {

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

    // get the read group and sample name
    uint8_t *rgptr = bam_aux_get(b, "RG");
    char* rg = (char*) (rgptr+1);
    //if (!rg_sample
    string sname;
    if (!rg_sample.empty()) {
        sname = rg_sample[string(rg)];
    }

    // add features to the alignment
    alignment.set_name(bam_get_qname(b));
    alignment.set_sequence(sequence);
    alignment.set_quality(quality);
    alignment.set_is_reverse(bam_is_rev(b));
    if (sname.size()) {
        alignment.set_sample_name(sname);
        alignment.set_read_group(rg);
    }

    return alignment;
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

void merge_alignments(Alignment& a1, const Alignment& a2) {
    int kept1, kept2;
    // determine likely order
    int64_t a1_start = a1.path().mapping(0).position().node_id();
    int64_t a2_start = a2.path().mapping(0).position().node_id();
    // and use this to change the direction of the merge
    // reverse?
    Path merged_path = a1_start <= a2_start ?
                                  merge_paths(a1.path(), a2.path(), kept1, kept2)
                                  : merge_paths(a2.path(), a1.path(), kept2, kept1);
    if (a1_start <= a2_start) {
        string seq = a1.sequence().substr(0, kept1) + a2.sequence().substr(a2.sequence().size()-kept2);
        a1.set_sequence(seq);
        *a1.mutable_path() = merged_path;
    } else {
        string seq = a2.sequence().substr(0, kept2) + a1.sequence().substr(a1.sequence().size()-kept1);
        a1.set_sequence(seq);
        *a1.mutable_path() = merged_path;
    }
    if (path_to_length(a1.path()) != a1.sequence().size()) {
        cerr << path_to_length(a1.path()) << " (to_length) != " << a1.sequence().size() << " (sequence size)" << endl;
        cerr << pb2json(a1) << endl;
        exit(1);
    }
}

void flip_nodes(Alignment& a, set<int64_t> ids, std::function<size_t(int64_t)> node_length) {
    Path* path = a.mutable_path();
    for(size_t i = 0; i < path->mapping_size(); i++) {
        // Grab each mapping
        Mapping* mapping = path->mutable_mapping(i);
        // Grab the position of each mapping
        Position* pos = mapping->mutable_position();
        
        if(ids.count(pos->node_id())) {
            // We need to flip this mapping
            
            // Flip its orientation
            mapping->set_is_reverse(!mapping->is_reverse());
            
            // Update the offset to count from the other end of the node.
            pos->set_offset(node_length(pos->node_id()) - pos->offset() - 1);
        } 
    }
    
}

}
