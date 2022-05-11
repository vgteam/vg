#ifndef VG_TEST_ALIGNER_HPP_INCLUDED
#define VG_TEST_ALIGNER_HPP_INCLUDED

#include "../aligner.hpp"

namespace vg {
namespace unittest {

// A silly shim to expose aligners for unittests now
// that I've made them inconvenient to use, except through
// the AlignerClient, which doesn't have public methods

/// We define a child class to expose all the protected stuff for testing
class TestAligner : public AlignerClient {
public:
    TestAligner(double gc_content) : AlignerClient(gc_content) {}
    TestAligner() : AlignerClient() {}
    using AlignerClient::set_alignment_scores;
    using AlignerClient::get_qual_adj_aligner;
    using AlignerClient::get_regular_aligner;
};

}
}

#endif
