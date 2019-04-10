//
//  packed_structs.hpp
//

#include "packed_structs.hpp"


namespace vg {
    
const double PackedVector::factor = 1.25;
const double PackedDeque::factor = 1.25;

PackedVector::PackedVector() {
    vec.width(1); // by default we start as a bitvector
}
    
PackedVector::PackedVector(istream& in) {
    deserialize(in);
}
    
PackedVector::~PackedVector() {
    
}
    
void PackedVector::deserialize(istream& in) {
    sdsl::read_member(filled, in);
    vec.load(in);
}
    
void PackedVector::serialize(ostream& out) const {
    sdsl::write_member(filled, out);
    vec.serialize(out);
}
    
PagedVector::PagedVector(size_t page_size) : page_size(page_size) {
    
}

PagedVector::PagedVector(istream& in) {
    deserialize(in);
}

PagedVector::~PagedVector() {
    
}
    
void PagedVector::deserialize(istream& in) {
    sdsl::read_member(filled, in);
    sdsl::read_member(page_size, in);
    anchors.deserialize(in);
    for (size_t i = 0; i < anchors.size(); i++) {
        pages.emplace_back(in);
    }
}

void PagedVector::serialize(ostream& out) const  {
    sdsl::write_member(filled, out);
    sdsl::write_member(page_size, out);
    anchors.serialize(out);
    for (size_t i = 0; i < pages.size(); i++) {
        pages[i].serialize(out);
    }
}

PackedDeque::PackedDeque() {
    
}
    
PackedDeque::PackedDeque(istream& in) {
    deserialize(in);
}

PackedDeque::~PackedDeque() {
    
}
    
void PackedDeque::deserialize(istream& in) {
    sdsl::read_member(begin_idx, in);
    sdsl::read_member(filled, in);
    vec.deserialize(in);
}

void PackedDeque::serialize(ostream& out) const  {
    sdsl::write_member(begin_idx, out);
    sdsl::write_member(filled, out);
    vec.serialize(out);
}

}
