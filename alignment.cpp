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

void parse_rg_sample_map(char* hts_header, map<string, string>& rg_sample) {
    string header(hts_header);
    vector<string> header_lines = split(header, '\n');

    for (auto& line : header_lines) {

        // get next line from header, skip if empty
        if ( line.empty() ) { continue; }

        // lines of the header look like:
        // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
        //                     ^^^^^^^\ is our sample name
        if (line.find("@RG") == 0) {
            vector<string> rg_parts = split(line, "\t ");
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
    const string& quality = alignment.quality();
    string squality;
    std::transform(quality.begin(), quality.end(), squality.begin(),
                   [](char c) { return (char)quality_short_to_char((int8_t)c); });
    alignment.set_quality(squality);
}

void alignment_quality_char_to_short(Alignment& alignment) {
    const string& quality = alignment.quality();
    string squality;
    std::transform(quality.begin(), quality.end(), squality.begin(),
                   [](char c) { return (int8_t)quality_char_to_short((char)c); });
    alignment.set_quality(squality);
}

Alignment bam_to_alignment(const bam1_t *b, map<string, string>& rg_sample) {

    Alignment alignment;

    // get the sequence and qual
    int32_t lqseq = b->core.l_qseq;
    string sequence; sequence.resize(lqseq);
    string quality; quality.resize(lqseq);

    // take the quality as-is
    // we have to copy this way--- it's not null terminated
    uint8_t* qualptr = bam_get_qual(b);
    memcpy((char*)quality.c_str(), qualptr, lqseq);

    // process the sequence into chars
    uint8_t* seqptr = bam_get_seq(b);
    for (int i = 0; i < lqseq; ++i) {
        sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }

    // get the read group and sample name
    uint8_t *rgptr = bam_aux_get(b, "RG");
    char* rg = (char*) (rgptr+1);
    string& sname = rg_sample[string(rg)];

    // add features to the alignment
    alignment.set_name(bam_get_qname(b));
    alignment.set_sequence(sequence);
    alignment.set_quality(quality);
    alignment.set_is_reverse(bam_is_rev(b));
    alignment.set_sample_name(sname);
    alignment.set_read_group(rg);

    return alignment;
}

bool project_alignment(Alignment& alignment, Index& index, string& path_name,
                       string& chrom, int64_t& pos, string& cigar) {
    // establish a position and a mapping for the alignment against the other path
    Alignment projected;
    bool result = index.project_path(alignment.path(), path_name, projected);
    if (result) {
        
    }
}

int to_length(Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

int from_length(Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;

}

}
