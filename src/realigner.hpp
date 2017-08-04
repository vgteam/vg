#ifndef VG_REALIGNER_HPP_INCLUDED
#define VG_REALIGNER_HPP_INCLUDED

#include <iostream>
#include <map>
#include "vg.hpp"
#include "mapper.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "json2pb.h"

namespace vg {

using namespace std;

class Realigner {

public:

    Realigner(vcflib::VariantCallFile& v,
              FastaReference& r,
              const string& t)
        : reference(r)
        , vcf_file(v)
        , target(t)
        , mapper(nullptr)
        , gcsaidx(nullptr)
        , lcpid(nullptr)
        , xgidx(nullptr)
        , identity_trigger(0.9)
        , realign_unpaired(true)
        , softclip_trigger(2)
        , debug(false)
        , idx_kmer_size(16)
        , doubling_steps(3)
        , edge_max(0)
        , idx_path_only(false)
        {
            parse_region(target,
                         seq_name,
                         start_pos,
                         end_pos);
            construct();
        }

    FastaReference& ref;
    vcflib::VariantCallFile vcf_file;
    string target;
    string seq_name;
    int start_pos;
    int end_pos;

    bool debug;
    double identity_trigger;
    bool realign_unpaired;
    double softclip_trigger;

    int idx_kmer_size;
    int edge_max;
    bool idx_path_only;
    int doubling_steps;

    /*
    // todo
    int context_depth;
    int thread_extension;
    int max_attempts;
    double min_identity;
    int alignment_threads;
    int max_mem_length;
    int min_mem_length;
    int hit_max;
    bool max_target_factor;
    int max_multimaps;
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    */

    Mapper* mapper;
    gcsa::GCSA* gcsaidx;
    gcsa::LCPArray* lcpidx;
    xg::XG* xgidx;

    void construct(void);
    Alignment realign(const Alignment& aln);

};

}

#endif
