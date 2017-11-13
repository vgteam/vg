#include "counter.hpp"

namespace vg {

Counter::Counter(void) : xgidx(nullptr) { }
Counter::Counter(xg::XG* xidx) : xgidx(xidx) { }
Counter::~Counter(void) { }

void Counter::load(const vector<string>& file_names) {
}

void Counter::write(const string& file_name) {
}

void Counter::compact(void) {
}

void Counter::add(const Alignment& aln) {
}

}
