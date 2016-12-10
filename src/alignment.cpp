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
        //cerr << string_quality_short_to_char(quality) << endl;
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
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
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
            bool got_anything = false;
#pragma omp critical (fastq_input)
            {
                if (more_data) {
                    got_anything = more_data = get_next_alignment_from_fastq(fp, buf, len, aln);
                    nLines++;
                }
            }
            if (got_anything) {
                lambda(aln);
            }
        }
    }
    for (auto& buf : bufs) {
        delete[] buf;
    }
    gzclose(fp);
    return nLines;
}

size_t fastq_paired_interleaved_for_each_parallel(string& filename, function<void(Alignment&, Alignment&)> lambda) {
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
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
            bool got_anything = false;
#pragma omp critical (fastq_input)
            {
                if (more_data) {
                    got_anything = more_data = get_next_interleaved_alignment_pair_from_fastq(fp, buf, len, mate1, mate2);
                    nLines++;
                }
            }
            if (got_anything) {
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
    if (!fp1) {
        cerr << "[vg::alignment.cpp] couldn't open " << file1 << endl; exit(1);
    }
    gzFile fp2 = (file2 != "-") ? gzopen(file2.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp2) {
        cerr << "[vg::alignment.cpp] couldn't open " << file2 << endl; exit(1);
    }
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
            bool got_anything = false;
#pragma omp critical (fastq_input)
            {
                if (more_data) {
                    got_anything = more_data = get_next_alignment_pair_from_fastqs(fp1, fp2, buf, len, mate1, mate2);
                    nLines++;
                }
            }
            if (got_anything) {
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
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
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
    delete buffer;
    return nLines;
}

size_t fastq_paired_two_files_for_each(string& file1, string& file2, function<void(Alignment&, Alignment&)> lambda) {
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
                         const bool refrev,
                         const string& cigar,
                         const string& mateseq,
                         const int32_t matepos,
                         const int32_t tlen) {

    assert(!sam_header.empty());
    string sam_file = "data:" + sam_header + alignment_to_sam(alignment, refseq, refpos, refrev, cigar, mateseq, matepos, tlen);
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
                        const bool refrev,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen) {
                        
    // Determine flags, using orientation.
    int32_t flags = sam_flag(alignment, refrev);
    
    stringstream sam;

    sam << (!alignment.name().empty() ? alignment.name() : "*") << "\t"
        << flags << "\t"
        << (refseq.empty() ? "*" : refseq) << "\t"
        << refpos + 1 << "\t"
        //<< (alignment.path().mapping_size() ? refpos + 1 : 0) << "\t" // positions are 1-based in SAM, 0 means unmapped
        << alignment.mapping_quality() << "\t"
        << (alignment.has_path() && alignment.path().mapping_size() ? cigar : "*") << "\t"
        << (mateseq == refseq ? "=" : mateseq) << "\t"
        << matepos + 1 << "\t"
        << tlen << "\t"
        // Make sure sequence always comes out in reference forward orientation by looking at the flags.
        << (!alignment.sequence().empty() ? (refrev ? reverse_complement(alignment.sequence()) : alignment.sequence()) : "*") << "\t";
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

string cigar_string(vector<pair<int, char> >& cigar) {
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

void mapping_cigar(const Mapping& mapping, vector<pair<int, char> >& cigar) {
    for (const auto& edit : mapping.edit()) {
        if (edit.from_length() && edit.from_length() == edit.to_length()) {
// *matches* from_length == to_length, or from_length > 0 and offset unset
            // match state
            cigar.push_back(make_pair(edit.from_length(), 'M'));
        } else {
            // mismatch/sub state
// *snps* from_length == to_length; sequence = alt
            if (edit.from_length() == edit.to_length()) {
                cigar.push_back(make_pair(edit.from_length(), 'M'));
            } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
                int32_t del = edit.from_length() - edit.to_length();
                int32_t eq = edit.to_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(del, 'D'));
            } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
                int32_t ins = edit.to_length() - edit.from_length();
                int32_t eq = edit.from_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(ins, 'I'));
            }
        }
    }
}

// act like the path this is against is the reference
// and generate an equivalent cigar
// Produces CIGAR in forward strand space of the reference sequence.
string cigar_against_path(const Alignment& alignment, bool on_reverse_strand) {
    vector<pair<int, char> > cigar;
    if (!alignment.has_path()) return "";
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
    if (cigar.front().second == 'I') {
        cigar.front().second = 'S';
    }
    if (cigar.back().second == 'I') {
        cigar.back().second = 'S';
    }
    
    return cigar_string(cigar);
}

int32_t sam_flag(const Alignment& alignment, bool on_reverse_strand) {
    int16_t flag = 0;

    if (alignment.score() == 0) {
        // unmapped
        flag |= BAM_FUNMAP;
    } else {
        // correctly aligned
        flag |= BAM_FPROPER_PAIR;
    }
    if (on_reverse_strand) {
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
    res.set_name(aln.name());
    res.set_score(aln.score());
    //cerr << "drop " << drop << " from start" << endl;
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

// merge that properly handles long indels
// assumes that alignments should line up end-to-end
Alignment merge_alignments(const vector<Alignment>& alns, bool debug) {

    if (alns.size() == 0) {
        Alignment aln;
        return aln;
    } else if (alns.size() == 1) {
        return alns.front();
    }

    // execute a serial merge
    // buliding up the alignment
    Alignment merged;

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
            extend_alignment(merged, aln, debug);
        }
    }
    *merged.mutable_path() = simplify(merged.path());
    return merged;
}


Alignment& extend_alignment(Alignment& a1, const Alignment& a2, bool debug) {
    //if (debug) cerr << "extending alignment " << endl << pb2json(a1) << endl << pb2json(a2) << endl;
    a1.set_sequence(a1.sequence() + a2.sequence());
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
    a3.set_sequence(a1.sequence() + a2.sequence());
    *a3.mutable_path() = concat_paths(a1.path(), a2.path());
    if (debug) cerr << "merged alignments, result is " << endl << pb2json(a3) << endl;
    return a3;
}

void translate_nodes(Alignment& a, const map<id_t, pair<id_t, bool> >& ids, const std::function<size_t(int64_t)>& node_length) {
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

Alignment simplify(const Alignment& a) {
    auto aln = a;
    *aln.mutable_path() = simplify(aln.path());
    return aln;
}

void write_alignment_to_file(const Alignment& aln, const string& filename) {
    ofstream out(filename);
    vector<Alignment> alnz = { aln };
    stream::write_buffered(out, alnz, 1);
    out.close();
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

}
