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
       << endl;
    return os;
}

}
