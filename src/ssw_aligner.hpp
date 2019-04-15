#ifndef VG_SSW_ALIGNER_HPP_INCLUDED
#define VG_SSW_ALIGNER_HPP_INCLUDED

#include <vector>
#include <set>
#include <string>
#include "ssw_cpp.h"
#include <vg/vg.pb.h>
#include "path.hpp"

namespace vg {


class SSWAligner {
public:

    SSWAligner(
        uint8_t _match = 1,
        uint8_t _mismatch = 4,
        uint8_t _gap_open = 6,
        uint8_t _gap_extension = 1)
        : match(_match)
        , mismatch(_mismatch)
        , gap_open(_gap_open)
        , gap_extension(_gap_extension) { }

    ~SSWAligner(void) { }

    uint8_t match;
    uint8_t mismatch;
    uint8_t gap_open;
    uint8_t gap_extension;

    // alignment functions
    Alignment align(const string& query, const string& ref);
    Alignment ssw_to_vg(const StripedSmithWaterman::Alignment& ssw_aln,
                        const string& query, const string& ref);
    void PrintAlignment(const StripedSmithWaterman::Alignment& alignment);

};

} // end namespace vg

#endif
