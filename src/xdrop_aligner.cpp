/**
 * \file xdrop_aliigner.cpp: contains implementation of XdropAligner
 */

#include "dozeu_interface.hpp"

// Configure dozeu:
// We want the full length bonus included
#ifndef DZ_FULL_LENGTH_BONUS
#define DZ_FULL_LENGTH_BONUS
#endif
// We want the non-quality adjusted versions of functions
#ifdef DZ_QUAL_ADJ
#undef DZ_QUAL_ADJ
#endif
// We require these particular values for this enum because we index arrays with it.
enum { MISMATCH = 1, MATCH = 2, INS = 3, DEL = 4 };
// Set dozeu's CIGAR codes to match our enum
#ifndef DZ_CIGAR_OP
#define DZ_CIGAR_OP 0x04030201
#endif

// To turn on debugging:
//#define DEBUG
//#define DZ_PRINT_VECTOR
 
#include <dozeu/dozeu.h>

using namespace vg;

XdropAligner::XdropAligner(const XdropAligner& other)
{
    *this = other;
}

XdropAligner& XdropAligner::operator=(const XdropAligner& other)
{
	if (this != &other) {

        if (dz) {
            dz_destroy(dz);
        }
        dz = dz_init(other.dz->matrix,
                     *((const uint16_t*) &other.dz->giv),
                     *((const uint16_t*) &other.dz->gev));
    }

	return *this;
}

XdropAligner::XdropAligner(XdropAligner&& other)
{
    *this = other;
}

XdropAligner& XdropAligner::operator=(XdropAligner&& other)
{
	if (this != &other) {
        if (dz) {
            dz_destroy(dz);
        }
        dz = other.dz;
        other.dz = nullptr;
    }

	return *this;
}

XdropAligner::XdropAligner(const int8_t* _score_matrix, int8_t _gap_open, int8_t _gap_extension)
{
    // xdrop aligner uses the parameterization where both gap open and gap extend
    // are added when opening a gap
    assert(_gap_open - _gap_extension >= 0);
    assert(_gap_extension > 0);
    dz = dz_init(_score_matrix, _gap_open - _gap_extension, _gap_extension);
}

XdropAligner::~XdropAligner(void)
{
    dz_destroy(dz);
}

dz_query_s* XdropAligner::pack_query_forward(const char* seq, const uint8_t* qual,
                                             int8_t full_length_bonus, size_t len) {
    return dz_pack_query_forward(dz, seq, full_length_bonus, len);
}

dz_query_s* XdropAligner::pack_query_reverse(const char* seq, const uint8_t* qual,
                                             int8_t full_length_bonus, size_t len) {
    return dz_pack_query_reverse(dz, seq, full_length_bonus, len);
}

const dz_forefront_s* XdropAligner::scan(const dz_query_s* query, const dz_forefront_s** forefronts,
                                         size_t n_forefronts, const char* ref, int32_t rlen,
                                         uint32_t rid, uint16_t xt) {
    return dz_scan(dz, query, forefronts, n_forefronts, ref, rlen, rid, xt);
}

const dz_forefront_s* XdropAligner::extend(const dz_query_s* query, const dz_forefront_s** forefronts,
                                           size_t n_forefronts, const char* ref, int32_t rlen,
                                           uint32_t rid, uint16_t xt) {
    return dz_extend(dz, query, forefronts, n_forefronts, ref, rlen, rid, xt);
}

dz_alignment_s* XdropAligner::trace(const dz_forefront_s* forefront) {
    return dz_trace(dz, forefront);
}

void XdropAligner::flush() {
    dz_flush(dz);
}

/**
 * end of xdrop_aligner.cpp
 */
