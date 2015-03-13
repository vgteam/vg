#include "alignment.hpp"

namespace vg {

int sam_for_each(string& filename, function<void(const Alignment&)>& lambda) {

    hFILE *fp = hopen(filename.c_str(), "r");
    samFile *in = hts_hopen(fp, filename.c_str(), "r");
    if (in == NULL) return 0;
    bam_hdr_t *hdr = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    while (sam_read1(in, hdr, b) >= 0) {
        lambda(bam_to_alignment(b));
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(in);
    if (fp && hclose(fp) < 0) {
        cerr << "htsfile: closing " << filename << " failed" << endl;
    }
    return 1;

}

Alignment bam_to_alignment(bam1_t* b) {

    Alignment alignment;

    uint32_t lqname = b->core.l_qname;
    string name; name.resize(lqname);
    int32_t lqseq = b->core.l_qseq;
    string sequence; sequence.resize(lqseq);
    string quality; quality.resize(lqseq);

    char* nameptr = bam_get_qname(b);
    memcpy((char*)name.c_str(), nameptr, lqname);

    uint8_t* qualptr = bam_get_qual(b);
    memcpy((char*)quality.c_str(), qualptr, lqseq);

    uint8_t* seqptr = bam_get_seq(b);
    for (int i = 0; i < lqseq; ++i) {
        //int8_t c = seqptr[i];
        const char singleBase = BAM_DNA_LOOKUP[ ( (seqptr[(i/2)] >> (4*(1-(i%2)))) & 0xf ) ];
        sequence[i] = singleBase;
    }

    alignment.set_name(name);
    alignment.set_sequence(sequence);
    alignment.set_quality(quality);
    return alignment;
}

}
