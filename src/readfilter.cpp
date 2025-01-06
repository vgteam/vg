#include "readfilter.hpp"

namespace vg {

using namespace std;
using namespace vg::io;

ostream& operator<<(ostream& os, const Counts& counts) {
    os << "Total Filtered:                " << counts.counts[Counts::FilterName::filtered] << " / "
       << counts.counts[Counts::FilterName::read] << endl
       << "Read Name Filter:              " << counts.counts[Counts::FilterName::wrong_name] << endl
       << "Subsequence Filter:            " << counts.counts[Counts::FilterName::subsequence] << endl
       << "Proper Pair Filter:            " << counts.counts[Counts::FilterName::proper_pair] << endl
       << "Unmapped Filter:               " << counts.counts[Counts::FilterName::unmapped] << endl
       << "refpos Contig Filter:          " << counts.counts[Counts::FilterName::wrong_refpos] << endl
       << "Feature Filter:                " << counts.counts[Counts::FilterName::excluded_feature] << endl
       << "Min Identity Filter:           " << counts.counts[Counts::FilterName::min_score] << endl
       << "Min Secondary Identity Filter: " << counts.counts[Counts::FilterName::min_sec_score] << endl
       << "Max Overhang Filter:           " << counts.counts[Counts::FilterName::max_overhang] << endl
       << "Min End Match Filter:          " << counts.counts[Counts::FilterName::min_end_matches] << endl
       << "Split Read Filter:             " << counts.counts[Counts::FilterName::split] << endl
       << "Repeat Ends Filter:            " << counts.counts[Counts::FilterName::repeat] << endl
       << "All Defrayed Filter:           " << counts.counts[Counts::FilterName::defray_all] << endl
       << "Min Quality Filter:            " << counts.counts[Counts::FilterName::min_mapq] << endl
       << "Min Base Quality Filter:       " << counts.counts[Counts::FilterName::min_base_qual] << endl
       << "Random Filter:                 " << counts.counts[Counts::FilterName::random] << endl
       << "Annotation Filter:             " << counts.counts[Counts::FilterName::annotation] << endl
       << "Incorrectly Mapped Filter:     " << counts.counts[Counts::FilterName::incorrectly_mapped] << endl
       << "Max Reads Filter:              " << counts.counts[Counts::FilterName::max_reads] << endl
       << endl;
    return os;
}

// for some reason these and only these methods become duplicate symbols if they're included
// in the header? i don't really understand why
template<>
bool ReadFilter<Alignment>::is_mapped(const Alignment& alignment) const {
    return alignment.path().mapping_size() != 0;
}

template<>
bool ReadFilter<MultipathAlignment>::is_mapped(const MultipathAlignment& mp_alignment) const {
    for (size_t i = 0; i < mp_alignment.subpath_size(); ++i) {
        if (mp_alignment.subpath(i).path().mapping_size() != 0) {
            return true;
        }
    }
    return false;
}


}
