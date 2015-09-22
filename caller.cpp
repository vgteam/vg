#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "caller.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

const double Caller::Log_zero = (double)-1e100;

// these values pretty arbitrary at this point
const int Caller::Default_buffer_size = 1000;
const double Caller::Default_het_prior = 0.001; // from MAQ
const int Caller::Default_min_depth = 5;
const double Caller::Default_min_frac = 0.75;
const double Caller::Default_min_likelihood = 1e-50;
const int Caller::Default_context = 5;
const char Caller::Default_default_quality = 10;

void Caller::clear() {
    _alignments.clear();
}

void Caller::to_json(ostream& out) {
    for (auto aln : _alignments) {
        out << pb2json(aln);
    }
}

void Caller::write(ostream& out) {
    function<Alignment&(uint64_t)> lambda = [&](uint64_t i) -> Alignment& {
        return _alignments[i];
    };
    stream::write(out, _alignments.size(), lambda);
}

void Caller::flush_buffer(ostream& out, bool json) {    
#pragma omp critical (out)
    {
        if (json) {
            to_json(out);
        } else {
            write(out);
        }
    }
    _alignments.clear();
}

void Caller::call_node_pileup(const NodePileup& pileup, ostream& out, bool json) {
    for (int i = 0; i < pileup.base_pileup_size(); ++i) {
        if (pileup.base_pileup(i).num_bases() >= _min_depth) {
            call_base_pileup(pileup, i);
        }
    }
    if (_alignments.size() >= _buffer_size) {
        flush_buffer(out, json);
    }
}

void Caller::call_base_pileup(const NodePileup& np, int64_t offset) {
    const BasePileup& bp = np.base_pileup(offset);
    
    // parse the pilueup structure
    vector<pair<int, int> > base_offsets;
    vector<pair<int, int> > indel_offsets;
    Pileups::parse_base_offsets(bp, base_offsets, indel_offsets, indel_offsets);

    // compute top two most frequent bases and their counts
    char top_base;
    int top_count;
    char second_base;
    int second_count;
    compute_top_frequencies(bp, base_offsets, top_base, top_count, second_base, second_count);

    // note first and second base will be upper case too
    char ref_base = ::toupper(bp.ref_base());

    // test against thresholding heuristics
    if (top_count + second_count >= _min_depth &&
        (double)(top_count + second_count) / (double)base_offsets.size() >= _min_frac) {

        // compute max likelihood snp genotype.  it will be one of the three combinations
        // of the top two bases (we don't care about case here)
        pair<char, char> g = mp_snp_genotype(bp, base_offsets, top_base, second_base);
        
        // if any of these bases are different than the reference, add them as snps
        if (g.first != ref_base) {
            create_snp(np, offset, g.first);
        }
        if (g.second != ref_base && g.second != g.first) {
            create_snp(np, offset, g.second);
        }
    }
}

void Caller::compute_top_frequencies(const BasePileup& bp,
                                     const vector<pair<int, int> >& base_offsets,
                                     char& top_base, int& top_count,
                                     char& second_base, int& second_count) {

    const string& bases = bp.bases();
    // frequency of each base (nidx / idxn converts to from / int)
    int hist[5] = {0};
    for (auto i : base_offsets) {
        char base = Pileups::extract_match(bp, i.first);
        ++hist[nidx(base)];
    }
    
    int first = max_element(hist, hist + 4) - hist;
    int second = first == 0 ? 1 : 0;
    for (int i = 0; i < 4; ++i) {
        if (i != first && hist[i] > hist[second]) {
            second = i;
        }
    }

    // break ties with reference
    int refidx = nidx(bp.ref_base());
    if (hist[first] == hist[refidx]) {
        first = refidx;
    } else if (hist[second] == hist[refidx]) {
        second = refidx;
    }

    top_base = idxn(first);
    top_count = hist[first];
    second_base = idxn(second);
    second_count = hist[second];
}

// Estimate the most probable snp genotype
pair<char, char> Caller::mp_snp_genotype(const BasePileup& bp,
                                         const vector<pair<int, int> >& base_offsets,
                                         char top_base, char second_base) {
    char ref_base = ::toupper(bp.ref_base());

    // gotta do better than this:
    pair<char, char> mp_genotype(ref_base, ref_base);
    double mp = _min_log_likelihood + _hom_log_prior;

    // genotype with 0 top_bases
    double gl = genotype_log_likelihood(bp, base_offsets, 0, top_base, second_base);
    double p = _hom_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(second_base, second_base);
    }

    // genotype with 1 top_base
    gl = genotype_log_likelihood(bp, base_offsets, 1, top_base, second_base);
    p = _het_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(top_base, second_base);
    }

    // genotype with 2 top_bases
    gl = genotype_log_likelihood(bp, base_offsets, 2, top_base, second_base);
    p = _hom_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(top_base, top_base);
    }
    
    // note, we're throwing away the probabilty value (mp) here.
    // should figure out where to stick it in the output. 
    return mp_genotype;
}

// This is Equation 2 (tranformed to log) from
// A statistical framework for SNP calling ... , Heng Li, Bioinformatics, 2011
// http://bioinformatics.oxfordjournals.org/content/27/21/2987.full
double Caller::genotype_log_likelihood(const BasePileup& bp,
                                       const vector<pair<int, int> >& base_offsets,
                                       double g, char first, char second) {
    double m = 2.; // always assume two alleles

    double log_likelihood = log(0.25); // 1 / m^2, where m = ploidy = 2;

    const string& bases = bp.bases();
    const string& quals = bp.qualities();
    double perr;

    for (int i = 0; i < base_offsets.size(); ++i) {
        char base = Pileups::extract_match(bp, base_offsets[i].first);
        char qual = base_offsets[i].second >= 0 ? quals[base_offsets[i].second] : _default_quality;
        perr = phred2prob(qual);
        if (base == first) {
            log_likelihood += safe_log((m - g) * perr + g * (1. - perr));
        } else if (base == second) {
            log_likelihood += safe_log((m - g) * (1. - perr) + g * perr);
        } else {
            log_likelihood += safe_log(perr * perr);
        }
    }

    return log_likelihood;
}


void Caller::create_snp(const NodePileup& np, int64_t offset, char base) {
    _alignments.push_back(Alignment());
    Alignment& alignment = _alignments.back();
    string sequence;
    // note on context: since we don't have access to the Node object,
    // need to get sequence out of the pileups...
    // left context (run right-to-left so we can cut at empty base more easily
    for (int64_t i = offset - 1; i >= std::max((int64_t)0, offset - _context); --i) {
        const BasePileup* bp = Pileups::get_base_pileup(np, i);
        if (!bp || bp->num_bases() == 0) {
            break;
        } else {
            sequence.insert(0, 1, (char)bp->ref_base());
        }
    }
    // snp base
    int seq_offset = sequence.length();
    sequence += base;
    // right context
    for (int64_t i = offset + 1;
         i < std::min((int64_t)np.base_pileup_size(), offset + 1 + _context); ++i) {
        const BasePileup* bp = Pileups::get_base_pileup(np, i);
        if (!bp || bp->num_bases() == 0) {
            break;
        } else {
            sequence += (char)bp->ref_base();
        }
    }
    alignment.set_sequence(sequence);
    Path* path = alignment.mutable_path();
    Mapping* mapping = path->add_mapping();
    mapping->mutable_position()->set_node_id(np.node_id());
    mapping->mutable_position()->set_offset(offset - seq_offset);
    mapping->set_is_reverse(false);
    // create match edit from left context
    if (seq_offset > 0) {
        Edit* edit = mapping->add_edit();
        edit->set_to_length(seq_offset);
        edit->set_from_length(seq_offset);
    }
    // create edit from snp
    Edit* edit = mapping->add_edit();
    edit->set_to_length(1);
    edit->set_from_length(1);
    edit->set_sequence(string(1, base));
    // create edit from right context
    if (seq_offset < sequence.length() - 1) {
        Edit* edit = mapping->add_edit();
        edit->set_to_length(sequence.length() - 1 - seq_offset);
        edit->set_from_length(sequence.length() - 1 - seq_offset);
    }
}


}
