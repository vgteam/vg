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

bool VcfBuffer::has_tabix() {
    return file && file->usingTabix;
}

bool VcfBuffer::set_region(const string& contig, int64_t start, int64_t end) {
    if(!has_tabix()) {
        // Nothing to seek in (or no file)
        return false;
    }

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

}
