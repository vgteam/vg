/**
 * \file
 * constructor.cpp: contains implementations for vg construction functions.
 */


#include "vcf_buffer.hpp"

namespace vg {

using namespace std;

vcflib::Variant* VcfBuffer::get() {
    if(has_buffer) {
        // We have a variant
        return &buffer;
    } else {
        // No next variant loaded.
        return nullptr;
    }
}

void VcfBuffer::handle_buffer() {
    // The variant in the buffer has been dealt with
    assert(has_buffer);
    has_buffer = false;
}

void VcfBuffer::fill_buffer() {
    if(file != nullptr && file->is_open() && !has_buffer && safe_to_get) {
        // Put a new variant in the buffer if we have a file and the buffer was empty.
        has_buffer = safe_to_get = file->getNextVariant(buffer);
        if(has_buffer) {
            // Convert to 0-based positions.
            // TODO: refactor to use vcflib zeroBasedPosition()...
            buffer.position -= 1;
        }
#ifdef debug
        cerr << "Variant in buffer: " << buffer << endl;
#endif
    }
}

bool VcfBuffer::has_tabix() const {
    return file && file->usingTabix;
}

bool VcfBuffer::set_region(const string& contig, int64_t start, int64_t end) {
    if(!has_tabix()) {
        // Nothing to seek in (or no file)
        return false;
    }

    // Discard any variants we had.
    has_buffer = false;
    
    // Remember that we can get the next variant now, in case we had hit the end
    // of the VCF.
    safe_to_get = true;

    if(start != -1 && end != -1) {
        // We have a start and end
        return file->setRegion(contig, start, end);
    } else {
        // Just seek to the whole chromosome
        return file->setRegion(contig);
    }
}

VcfBuffer::VcfBuffer(vcflib::VariantCallFile* file) : file(file) {
    // Our buffer needs to know about the VCF file it is reading from, because
    // it cares about the sample names. If it's not associated properely, we
    // can't getNextVariant into it.
    if (file) {
        // But only do it if we actually have a real file
        buffer.setVariantCallFile(file);
    }
}


WindowedVcfBuffer::WindowedVcfBuffer(vcflib::VariantCallFile* file, size_t window_size): reader(file), window_size(window_size) {
    // Nothing to do!
}

bool WindowedVcfBuffer::has_tabix() const {
    return reader.has_tabix();
}

bool WindowedVcfBuffer::set_region(const string& contig, int64_t start, int64_t end) {
    // Clear buffers
    variants_before.clear();
    variants_after.clear();
    has_current = false;
    
    return reader.set_region(contig, start, end);
}

bool WindowedVcfBuffer::next() {
    if (has_current) {
        // Put the current variant on the list of previous variants
        variants_before.emplace_back(std::move(current));
        has_current = false;
    }
    
    if (!variants_after.empty()) {
        // We already know what the new current variant should be, so grab it
        current = std::move(variants_after.front());
        variants_after.pop_front();
        has_current = true;
    } else {
        // We need to read the next variant from the VCF and make that the
        // current variant.
        reader.fill_buffer();
        if (reader.get() == nullptr) {
            // Couldn't find anything
            return false;
        } else {
            // Copy what we found
            current = *(reader.get());
            reader.handle_buffer();
            has_current = true;
        }
    }
    
    // We will always have a current variant if we get here.
    assert(has_current);
    
    // Now we know we have a current variant. We need to drop variants on the
    // left that are on the wrong sequence or have gotten too far away.
    while (!variants_before.empty() &&
        (current.sequenceName != variants_before.front().sequenceName ||
        current.position - variants_before.front().position > window_size)) {
        
        // As long as the leftmost variant exists and is on the wrong contig or
        // too far left, pop it.
        variants_before.pop_front();
    } 
    
    
    // We also need to load in new variants on the right until what's buffered
    // is too far away, or we run out.
    reader.fill_buffer();
    while (reader.get() != nullptr &&
        reader.get()->sequenceName == current.sequenceName &&
        reader.get()->position - current.position <= window_size) {
        
        // As long as we have a next variant on this contig that's in range,
        // grab a copy.
        variants_after.push_back(*(reader.get()));
        reader.handle_buffer();
        reader.fill_buffer();
    }
    
    // We have a (new) current variant
    return true;
}

tuple<vector<vcflib::Variant*>, vcflib::Variant*, vector<vcflib::Variant*>> WindowedVcfBuffer::get() {
    if (!has_current) {
        throw runtime_error("Can't get variant when no current variant is loaded.");
    }
    
    // We have to make vectors of all the pointers.
    vector<vcflib::Variant*> before;
    vector<vcflib::Variant*> after;

    transform(variants_before.begin(), variants_before.end(), back_inserter(before), [](vcflib::Variant& v) { return &v; });
    transform(variants_after.begin(), variants_after.end(), back_inserter(after), [](vcflib::Variant& v) { return &v; });
    
    return make_tuple(before, &current, after);
}

tuple<vector<vcflib::Variant*>, vcflib::Variant*, vector<vcflib::Variant*>> WindowedVcfBuffer::get_nonoverlapping() {
    if (!has_current) {
        throw runtime_error("Can't get variant when no current variant is loaded.");
    }
    
    // We have to make vectors of all the pointers.
    vector<vcflib::Variant*> before;
    vector<vcflib::Variant*> after;
    
    // What's the last base used by a variant? We need this to strip out
    // variants that overlap each other.
    int64_t last_position_used = 0;
    
    for (auto& v : variants_before) {
        // For every variant before the current one
        if (v.position > last_position_used && v.position + v.ref.size() <= current.position) {
            // The end of this variant is before the start of our current variant, so use it.
            before.push_back(&v);
            // Remember it uses its start base, plus 1 base per ref allele base
            last_position_used = v.position + v.ref.size() - 1;
        }
    }
    
    // Now mark used through the end of the current variant
    last_position_used = current.position + current.ref.size() - 1;
    
    for (auto& v : variants_after) {
        // For every variant after the current one
        if (v.position > last_position_used) {
            // The end of the last variant is before the start of this one, so use it.
            after.push_back(&v);
            // Remember it uses its start base, plus 1 base per ref allele base
            last_position_used = v.position + v.ref.size() - 1;
        }
    }
    
    return make_tuple(before, &current, after);
}

}





















