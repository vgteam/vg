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
        // keep trailing /1 /2
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
    if (alignment.is_secondary()) {
        flag |= BAM_FSECONDARY;
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
    alignment.set_sequence(sequence);
    alignment.set_quality(quality);
    alignment.set_is_reverse(bam_is_rev(b));
    // TODO: htslib doesn't wrap this flag for some reason.
    alignment.set_is_secondary(b->core.flag & BAM_FSECONDARY);
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

Alignment strip_from_start(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    cerr << "drop " << drop << " from start" << endl;
    res.set_sequence(aln.sequence().substr(drop));
    if (!aln.has_path()) return res;
    *res.mutable_path() = cut_path(aln.path(), drop).second;
    assert(res.has_path());
    if (alignment_to_length(res) != res.sequence().size()) {
        cerr << "failed!!! drop from start 轰" << endl;
        cerr << pb2json(res) << endl << endl;
        assert(false);
    }
    return res;
}

Alignment strip_from_end(const Alignment& aln, size_t drop) {
    if (!drop) return aln;
    Alignment res;
    cerr << "drop " << drop << " from end" << endl;
    size_t cut_at = aln.sequence().size()-drop;
    cerr << "Cut at " << cut_at << endl;
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

// merge that properly handles long indels
// assumes that alignments should line up end-to-end
Alignment merge_alignments(const vector<Alignment>& alns, const vector<size_t>& overlaps, bool debug) {
    if (alns.size() == 0) {
        Alignment aln;
        return aln;
    } else if (alns.size() == 1) {
        return alns.front();
    }
    Alignment aln = alns.front(); // keep the whole first alignment
    for (size_t i = 1; i < alns.size(); ++i) {
        aln = merge_alignments(aln, alns[i], overlaps[i], debug);
    }
    return aln;
}

Alignment merge_alignments(Alignment a1, Alignment a2, size_t overlap, bool debug) {
    //cerr << "overlap is " << overlap << endl;
    // if either doesn't have a path, then treat it like a massive softclip
    if (!a1.has_path()) {
        Mapping m;
        Edit* e = m.add_edit();
        e->set_to_length(a1.sequence().size());
        e->set_sequence(a1.sequence());
        *a1.mutable_path()->add_mapping() = m;
    }
    if (!a2.has_path()) {
        Mapping m;
        Edit* e = m.add_edit();
        e->set_to_length(a2.sequence().size());
        e->set_sequence(a2.sequence());
        *a2.mutable_path()->add_mapping() = m;
    }
    if (debug) cerr << "merging alignments " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    // splice em up in the middle
    a1 = strip_from_end(a1, overlap/2);
    a2 = strip_from_start(a2, overlap-(overlap/2));
    // and concatenate them
    Alignment a3;
    a3.set_sequence(a1.sequence() + a2.sequence());
    *a3.mutable_path() = concat_paths(a1.path(), a2.path());
    if (debug) cerr << "merged alignments " << pb2json(a3) << endl;
    return a3;    
}

// note only handles alignments against the forward strand
/*
Alignment wild_merge_alignments(Alignment a1, Alignment a2, size_t overlap, bool debug) {
    //cerr << "overlap is " << overlap << endl;
    // if either doesn't have a path, then treat it like a massive softclip
    if (!a1.has_path()) {
        Mapping m;
        Edit* e = m.add_edit();
        e->set_to_length(a1.sequence().size());
        e->set_sequence(a1.sequence());
        *a1.mutable_path()->add_mapping() = m;
        // todo we should check if a2 has a position
        // if so we'll add a1 on as a big soft clip
    }
    if (!a2.has_path()) {
        Mapping m;
        Edit* e = m.add_edit();
        e->set_to_length(a2.sequence().size());
        e->set_sequence(a2.sequence());
        *a2.mutable_path()->add_mapping() = m;
        // todo, check if a1 has a position
        // if so, add a2 on as a soft clip
        // gaaaaa
    }
    if (debug) cerr << "merging alignments " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    // keep things simple and assume that a1 is aligned "left" of a2
    // determine likely order
    auto& a1_first_mapping = a1.path().mapping(0);
    auto& a1_last_mapping =  a1.path().mapping(a1.path().mapping_size()-1);
    auto& a2_first_mapping = a2.path().mapping(0);
    auto& a2_last_mapping =  a2.path().mapping(a2.path().mapping_size()-1);
    Position a1_start_pos =  a1_first_mapping.position();
    Position a1_end_pos =    a1_last_mapping.position();
    // set position to point to 1-past end
    a1_end_pos.set_offset(from_length(a1_last_mapping));
    Position a2_start_pos =  a2_first_mapping.position();
    Position a2_end_pos =    a2_last_mapping.position();
    // set position to point to 1-past end
    a2_end_pos.set_offset(from_length(a2_last_mapping));

    // find a common  position
    Position common_pos;
    bool found_common_pos = false;
    
    // iterate forwards from start of a2 looking for a place that approximately
    // balances the overlap between both and is not in a soft clip or indel
    // ....
    for (size_t i = 0; i < a2.path().mapping_size(); ++i) {
        const Mapping& m = a2.path().mapping(i);
        Position p = m.position();
        // we can be more efficient by checking if we actually have mappings
        // for both paths against this node
        if (!maps_to_node(a1.path(), p.node_id())
            || !maps_to_position(a2.path(), p.node_id())) {
            continue;
        }
        // if they do, let's see if we can find an overlap between them
        for (size_t j = 0; j < mapping_from_length(m); ++j) {
            auto a1_splits = cut_path(a1.path(), m.position());
            auto& a1l = a1_splits.first;
            auto& a1r = a1_splits.second;
            auto a2_splits = cut_path(a2.path(), m.position());
            auto& a2l = a2_splits.first;
            auto& a2r = a2_splits.second;
            // check if the left path ends before this position
            // and break if so -- this suggests a gap and we can't merge
            auto& a1l_last = a1l.mapping(a1l.mapping_size()-1);
            if (a1l_last.position().offset() + from_length(a1l_last) < j) break;
            // check if the right path starts after this position
            // and continue if so
            auto& a2r_first = a2r.mapping(0);
            if (a2r_first.position().offset() > j) continue;
            // if we make it here it means we're overlapping or abutting at this point
            // although we may not meet other requirements, so verify before calling
            // this position + j our merge point

            // does the left side of the last path end in a softclip here?
            // does the right side of the next path start in a softclip?

            // do we meet our objectives for the overlap of the two alignments?
            
            
            // not before
            //size_t a1_after = to_length_after_pos(a1, m.position());
            //size_t a2_before = to_length_before_pos(a2, m.position());
        }
        // if we're in a match, at the start of a node, and == to the other mapping
        //cerr << "mapping " << pb2json(m) << endl;
        cerr << "evaluating possible cut position " << pb2json(m.position()) << endl;;
        // show the cuts
        // this shit ain't working

        
        cerr << "a1_after " << a1_after << " a2_before " << a2_before << endl;
        if (a1_after > 0 && a2_before > 0
            && a2_before >= overlap/4
            //&& m.position().offset() == 0
            && edit_is_match(m.edit(0))) {
            common_pos = m.position();
            found_common_pos = true;
            break;
        }
    }

    // remove the overlap from the alignments
    // so that we could concatenate them to get a new valid alignment
    if (!found_common_pos) {
        cerr << "could not find common position!" << endl;
        cerr << "must be a big gap" << endl;
        cerr << "to strip " << overlap/2 << " from end of first and " << overlap-(overlap/2) << " from second" << endl;
        a1 = strip_from_end(a1, overlap/2);
        a2 = strip_from_start(a2, overlap-(overlap/2));
        cerr << "first after strip " << pb2json(a1) << endl;
        cerr << "second after strip " << pb2json(a2) << endl;
    } else {
        cerr << "found a common position " << pb2json(common_pos) << endl;
        int a1_drop_to_common = to_length_after_pos(a1, common_pos);
        int a2_drop_to_common = to_length_before_pos(a2, common_pos);
        // count bp before our common position in each alignment
        a1 = strip_from_end(a1, a1_drop_to_common);
        a2 = strip_from_start(a2, a2_drop_to_common);
    }

    // how many bp of soft clips do we have on each side of our merge
    int left_softclip = softclip_end(a1);
    int right_softclip = softclip_start(a2);

    // possibilities
    
    // we are overlapping and there are no soft clips
    // ---- suggests clean mapping, easy to merge, remove overlap/2 from each bit
    
    // we are overlapping but there are soft clips
    // ---- prefer removal of soft clipping bits, then trim down evenly
    //      gap after this may indicate a smaller del
    
    // we are not overlapping and there are soft clips
    // ---- sounds like a big deletion. trim off the overlap and merge.
    
    // we are not overlapping and there are not soft clips (weird)
    // ---- this one is strange; it suggests some kind of deletion event or mismapping
    //      we should send a warning to the console, trim the overlap and be on our way
    
    Alignment a3;
    a3.set_sequence(a1.sequence() + a2.sequence());
    *a3.mutable_path() = concat_paths(a1.path(), a2.path());
    if (debug) cerr << "merged alignments " << pb2json(a3) << endl;
    return a3;
}
*/

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

int softclip_start(Alignment& alignment) {
    if (alignment.mutable_path()->mapping_size() > 0) {
        Path* path = alignment.mutable_path();
        Mapping* first_mapping = path->mutable_mapping(0);
        Edit* first_edit = first_mapping->mutable_edit(0);
        if (first_edit->from_length() == 0 && first_edit->to_length() > 0) {
            return first_edit->to_length();
        }
    }
    return 0;
}

int softclip_end(Alignment& alignment) {
    if (alignment.mutable_path()->mapping_size() > 0) {
        Path* path = alignment.mutable_path();
        Mapping* last_mapping = path->mutable_mapping(path->mapping_size()-1);
        Edit* last_edit = last_mapping->mutable_edit(last_mapping->edit_size()-1);
        if (last_edit->from_length() == 0 && last_edit->to_length() > 0) {
            return last_edit->to_length();
        }
    }
    return 0;
}

size_t to_length_after_pos(const Alignment& aln, const Position& pos) {
    cerr << "Getting to length after " << pb2json(pos) << endl;
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


}
