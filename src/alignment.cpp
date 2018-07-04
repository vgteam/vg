#include "alignment.hpp"
#include "stream.hpp"

#include <regex>

namespace vg {

int hts_for_each(string& filename, function<void(Alignment&)> lambda, xg::XG* xgindex) {

    samFile *in = hts_open(filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    map<string, string> rg_sample;
    parse_rg_sample_map(hdr->text, rg_sample);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        Alignment a = bam_to_alignment(b, rg_sample, hdr, xgindex);
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

int hts_for_each_parallel(string& filename, function<void(Alignment&)> lambda, xg::XG* xgindex) {

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
                Alignment a = bam_to_alignment(b, rg_sample, hdr, xgindex);
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
    if (0!=gzgets(fp,buffer,len)) {
        buffer[strlen(buffer)-1] = '\0';
        string name = buffer;
        if (name[0] == '@') {
            is_fasta = false;
        } else if (name[0] = '>') {
            is_fasta = true;
        } else {
            throw runtime_error("Found unexpected delimiter " + name.substr(0,1) + " in fastq/fasta input");
        }
        name = name.substr(1, name.find(' ')); // trim off leading @ and things after the first whitespace
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
    if (!is_fasta) {
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
    }

    return true;

}

bool get_next_interleaved_alignment_pair_from_fastq(gzFile fp, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp, buffer, len, mate1) && get_next_alignment_from_fastq(fp, buffer, len, mate2);
}

bool get_next_alignment_pair_from_fastqs(gzFile fp1, gzFile fp2, char* buffer, size_t len, Alignment& mate1, Alignment& mate2) {
    return get_next_alignment_from_fastq(fp1, buffer, len, mate1) && get_next_alignment_from_fastq(fp2, buffer, len, mate2);
}

size_t unpaired_for_each_parallel(function<bool(Alignment&)> get_read_if_available, function<void(Alignment&)> lambda) {

    size_t nLines = 0;
    vector<Alignment> *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_read_if_available, lambda)
#pragma omp single
    {
        
        // number of reads in each batch
        const uint64_t batch_size = 1 << 9; // 512
        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = 1 << 9; // 512
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
                                           function<bool(void)> single_threaded_until_true) {
    
    
    size_t nLines = 0;
    vector<pair<Alignment, Alignment> > *batch = nullptr;
    // number of batches currently being processed
    uint64_t batches_outstanding = 0;
    
#pragma omp parallel default(none) shared(batches_outstanding, batch, nLines, get_pair_if_available, single_threaded_until_true, lambda)
#pragma omp single
    {
        
        // number of pairs in each batch
        const uint64_t batch_size = 1 << 9; // 512
        // max # of such batches to be holding in memory
        uint64_t max_batches_outstanding = 1 << 9; // 512
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

size_t fastq_unpaired_for_each_parallel(const string& filename, function<void(Alignment&)> lambda) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 2 << 22; // 4M
    char* buf = new char[len];
    
    function<bool(Alignment&)> get_read = [&](Alignment& aln) {
        return get_next_alignment_from_fastq(fp, buf, len, aln);;
    };
    
    
    size_t nLines = unpaired_for_each_parallel(get_read, lambda);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
    
}

size_t fastq_paired_interleaved_for_each_parallel(const string& filename, function<void(Alignment&, Alignment&)> lambda) {
    return fastq_paired_interleaved_for_each_parallel_after_wait(filename, lambda, [](void) {return true;});
}
    
size_t fastq_paired_two_files_for_each_parallel(const string& file1, const string& file2, function<void(Alignment&, Alignment&)> lambda) {
    return fastq_paired_two_files_for_each_parallel_after_wait(file1, file2, lambda, [](void) {return true;});
}
    
size_t fastq_paired_interleaved_for_each_parallel_after_wait(const string& filename,
                                                             function<void(Alignment&, Alignment&)> lambda,
                                                             function<bool(void)> single_threaded_until_true) {
    
    gzFile fp = (filename != "-") ? gzopen(filename.c_str(), "r") : gzdopen(fileno(stdin), "r");
    if (!fp) {
        cerr << "[vg::alignment.cpp] couldn't open " << filename << endl; exit(1);
    }
    
    size_t len = 1 << 18; // 256k
    char* buf = new char[len];
    
    function<bool(Alignment&, Alignment&)> get_pair = [&](Alignment& mate1, Alignment& mate2) {
        return get_next_interleaved_alignment_pair_from_fastq(fp, buf, len, mate1, mate2);
    };
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true);
    
    delete[] buf;
    gzclose(fp);
    return nLines;
}
    
size_t fastq_paired_two_files_for_each_parallel_after_wait(const string& file1, const string& file2,
                                                           function<void(Alignment&, Alignment&)> lambda,
                                                           function<bool(void)> single_threaded_until_true) {
    
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
    
    size_t nLines = paired_for_each_parallel_after_wait(get_pair, lambda, single_threaded_until_true);
    
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

// Internal conversion function for both paired and unpaired codepaths
string alignment_to_sam_internal(const Alignment& alignment,
                                 const string& refseq,
                                 const int32_t refpos,
                                 const bool refrev,
                                 const string& cigar,
                                 const string& mateseq,
                                 const int32_t matepos,
                                 const int32_t tlen,
                                 bool paired) {
                        
    // Determine flags, using orientation, next/prev fragments, and pairing status.
    int32_t flags = sam_flag(alignment, refrev, paired);
   
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
    
    // We've observed some reads with the unmapped flag set and also a CIGAR string set, which shouldn't happen.
    // We will check for this. The CIGAR string will only be set in the output if the alignment has a path.
    assert((bool)(flags & BAM_FUNMAP) != (alignment.has_path() && alignment.path().mapping_size()));
    
    stringstream sam;
    
    string alignment_name;
    if (paired) {
        // We need to strip the /1 and /2 or _1 and _2 from paired reads so the two ends have the same name.
        alignment_name = regex_replace(alignment.name(), regex("[/_][12]$"), "");
    } else {
        // Keep the alignment name as is because even if the name looks paired, the reads are semantically unpaired.
        alignment_name = alignment.name();
    }
    
    sam << (!alignment_name.empty() ? alignment_name : "*") << "\t"
        << flags << "\t"
        << (mapped ? refseq : "*") << "\t"
        << refpos + 1 << "\t"
        << (mapped ? alignment.mapping_quality() : 0) << "\t"
        << (mapped ? cigar : "*") << "\t"
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
    sam << "\n";
    return sam.str();
}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, mateseq, matepos, tlen, true);

}

string alignment_to_sam(const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar) {
    
    return alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, "", -1, 0, false);

}

// Internal conversion function for both paired and unpaired codepaths
bam1_t* alignment_to_bam_internal(const string& sam_header,
                                  const Alignment& alignment,
                                  const string& refseq,
                                  const int32_t refpos,
                                  const bool refrev,
                                  const string& cigar,
                                  const string& mateseq,
                                  const int32_t matepos,
                                  const int32_t tlen,
                                  bool paired) {

    assert(!sam_header.empty());
    
    // Make a tiny SAM file. Remember to URL-encode it, since it may contain '%'
    string sam_file = "data:," + percent_url_encode(sam_header +
        alignment_to_sam_internal(alignment, refseq, refpos, refrev, cigar, mateseq, matepos, tlen, paired));
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

bam1_t* alignment_to_bam(const string& sam_header,
                        const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar,
                        const string& mateseq,
                        const int32_t matepos,
                        const int32_t tlen) {
    
    return alignment_to_bam_internal(sam_header, alignment, refseq, refpos, refrev, cigar, mateseq, matepos, tlen, true);

}

bam1_t* alignment_to_bam(const string& sam_header,
                        const Alignment& alignment,
                        const string& refseq,
                        const int32_t refpos,
                        const bool refrev,
                        const string& cigar) {
    
    return alignment_to_bam_internal(sam_header, alignment, refseq, refpos, refrev, cigar, "", -1, 0, false);

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
            //cerr << "match " << edit.from_length() << endl;
        } else {
            // mismatch/sub state
// *snps* from_length == to_length; sequence = alt
            if (edit.from_length() == edit.to_length()) {
                cigar.push_back(make_pair(edit.from_length(), 'M'));
                //cerr << "match " << edit.from_length() << endl;
            } else if (edit.from_length() > edit.to_length()) {
// *deletions* from_length > to_length; sequence may be unset or empty
                int32_t del = edit.from_length() - edit.to_length();
                int32_t eq = edit.to_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(del, 'D'));
                //cerr << "del " << edit.from_length() - edit.to_length() << endl;
            } else if (edit.from_length() < edit.to_length()) {
// *insertions* from_length < to_length; sequence contains relative insertion
                int32_t ins = edit.to_length() - edit.from_length();
                int32_t eq = edit.from_length();
                if (eq) cigar.push_back(make_pair(eq, 'M'));
                cigar.push_back(make_pair(ins, 'I'));
                //cerr << "ins " << edit.to_length() - edit.from_length() << endl;
            }
        }
    }
}

int64_t cigar_mapping(const bam1_t *b, Mapping* mapping, xg::XG* xgindex) {
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

void mapping_against_path(Alignment& alignment, const bam1_t *b, char* chr, xg::XG* xgindex, bool on_reverse_strand) {

    if (b->core.pos == -1) return;

    Mapping mapping;

    int64_t length = cigar_mapping(b, &mapping, xgindex);

    Alignment aln = xgindex->target_alignment(chr, b->core.pos, b->core.pos + length, "", on_reverse_strand, mapping);

    *alignment.mutable_path() = aln.path();

    Position* refpos = alignment.add_refpos();
    refpos->set_name(chr);
    refpos->set_offset(b->core.pos);
    refpos->set_is_reverse(on_reverse_strand);
}

// act like the path this is against is the reference
// and generate an equivalent cigar
// Produces CIGAR in forward strand space of the reference sequence.
string cigar_against_path(const Alignment& alignment, bool on_reverse_strand, int64_t& pos, size_t path_len, size_t softclip_suppress) {
    vector<pair<int, char> > cigar;
    if (!alignment.has_path() || alignment.path().mapping_size() == 0) return "";
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

    return cigar_string(cigar);
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
    } else if (flag & BAM_FPAIRED) {
        // Aligned and in a pair, so assume it's properly paired.
        // TODO: this relies on us not emitting improperly paired reads
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

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample, const bam_hdr_t *bh, xg::XG* xgindex) {

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
    
    if (xgindex != nullptr && bh != nullptr) {
        alignment.set_mapping_quality(b->core.qual);
        mapping_against_path(alignment, b, bh->target_name[b->core.tid], xgindex, b->core.flag & BAM_FREVERSE);
    }
    
    // TODO: htslib doesn't wrap this flag for some reason.
    alignment.set_is_secondary(b->core.flag & BAM_FSECONDARY);
    if (sname.size()) {
        alignment.set_sample_name(sname);
        alignment.set_read_group(rg);
    }

    return alignment;
}

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample) {
    return bam_to_alignment(b, rg_sample, nullptr, nullptr);
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
                       xg::XG* xgindex,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
    if (!bedstream) {
        cerr << "Unable to open bed file." << endl;
        return;
    }
    string row;
    string seq;
    size_t sbuf;
    size_t ebuf;
    string name;
    size_t score;
    string strand;

    for (int line = 1; getline(bedstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        ss >> seq;
        ss >> sbuf;
        ss >> ebuf;

        if (ss.fail()) {
            cerr << "Error parsing bed line " << line << ": " << row << endl;
        } else {
            ss >> name;
            assert(sbuf < ebuf);
            ss >> score;
            ss >> strand;

            bool is_reverse = false;
            if(!ss.fail() && strand.compare("-") == 0) {
                is_reverse = true;
            }

            if (xgindex->path_rank(seq) == 0) {
                // This path doesn't exist, and we'll get a segfault or worse if
                // we go look for positions in it.
                cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
            } else {
                Alignment alignment = xgindex->target_alignment(seq, sbuf, ebuf, name, is_reverse);

                out_alignments->push_back(alignment);
            }
        }
    }
}

void parse_gff_regions(istream& gffstream,
                       xg::XG* xgindex,
                       vector<Alignment>* out_alignments) {
    out_alignments->clear();
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

    for (int line = 1; getline(gffstream, row); ++line) {
        if (row.size() < 2 || row[0] == '#') {
            continue;
        }
        istringstream ss(row);
        getline(ss, seq, '\t');
        getline(ss, source, '\t');
        getline(ss, type, '\t');
        getline(ss, buf, '\t');
        sbuf = atoi(buf.c_str());
        getline(ss, buf, '\t');
        ebuf = atoi(buf.c_str());

        if (ss.fail() || !(sbuf < ebuf)) {
            cerr << "Error parsing gtf/gff line " << line << ": " << row << endl;
        } else {
            getline(ss, score, '\t');
            getline(ss, strand, '\t');
            getline(ss, num, '\t');
            getline(ss, annotations, '\t');
            vector<string> vals = split(annotations, ";");
            for (auto& s : vals) {
                if (s.find("Name=") == 0) {
                    name = s.substr(5);
                }
            }

            bool is_reverse = false;
            if(!ss.fail() && strand.compare("-") == 0) {
                is_reverse = true;
            }

            if (xgindex->path_rank(seq) == 0) {
                // This path doesn't exist, and we'll get a segfault or worse if
                // we go look for positions in it.
                cerr << "warning: path \"" << seq << "\" not found in index, skipping" << endl;
            } else {
                Alignment alignment = xgindex->target_alignment(seq, sbuf, ebuf, name, is_reverse);

                out_alignments->push_back(alignment);
            }
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

void alignment_set_distance_to_correct(Alignment& aln, const Alignment& base) {
    auto base_offsets = alignment_refpos_to_path_offsets(base);
    return alignment_set_distance_to_correct(aln, base_offsets);
}

void alignment_set_distance_to_correct(Alignment& aln, const map<string ,vector<pair<size_t, bool> > >& base_offsets) {
    auto aln_offsets = alignment_refpos_to_path_offsets(aln);
    // bail out if we can't compare
    if (!(aln_offsets.size() && base_offsets.size())) return;
    // otherwise find the minimum distance and relative orientation
    Position result;
    size_t min_distance = std::numeric_limits<size_t>::max();
    for (auto& path : aln_offsets) {
        auto& name = path.first;
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

}
