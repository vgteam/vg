#ifndef VG_VCF_BUFFER_HPP
#define VG_VCF_BUFFER_HPP

/** \file
 * vcf_buffer.hpp: defines a buffer for reading through VCF files with
 * look-ahead.
 */

// We need vcflib
#include "Variant.h"

namespace vg {

using namespace std;

/**
 * Provides a one-variant look-ahead buffer on a vcflib::VariantFile. Lets
 * construction functions peek and see if they want the next variant, or lets
 * them ignore it for the next construction function for a different contig to
 * handle. Ought not to be copied.
 *
 * Handles conversion from 1-based vcflib coordinates to 0-based vg coordinates.
 */
class VcfBuffer {

public:
    /**
     * Return a pointer to the buffered variant, or null if no variant is
     * buffered. Pointer is invalidated when the buffer is handled. The variant
     * will have a 0-based start coordinate.
     */
    vcflib::Variant* get();
    
    /**
     * To be called when the buffer is filled. Marks the buffered variant as
     * handled, discarding it, and allowing another to be read.
     */
    void handle_buffer();
    
    /**
     * Can be called when the buffer is filled or empty. If there is no variant
     * in the buffer, tries to load a variant into the buffer, if one can be
     * obtained from the file.
     */
    void fill_buffer();
    
    /**
     * This returns true if we have a tabix index, and false otherwise. If this
     * is false, set_region may be called, but will do nothing and return false.
     */
    bool has_tabix();
    
    /**
     * This tries to set the region on the underlying vcflib VariantCallFile to
     * the given contig and region, if specified. Coordinates coming in should
     * be 0-based,a nd will be converted to 1-based internally.
     *
     * Returns true if the region was successfully set, and false otherwise (for
     * example, if there is not tabix index, or if the given region is not part
     * of this VCF. Note that if there is a tabix index, and set_region returns
     * false, the position in the VCF file is undefined until the next
     * successful set_region call.
     *
     * If either of start and end are specified, then both of start and end must
     * be specified.
     */
    bool set_region(const string& contig, int64_t start = -1, int64_t end = -1);
    
    /**
     * Make a new VcfBuffer buffering the file at the given pointer (which must
     * outlive the buffer, but which may be null).
     */
    VcfBuffer(vcflib::VariantCallFile* file = nullptr);
    
protected:
    
    // This stores whether the buffer is populated with a valid variant or not
    bool has_buffer = false;
    // This stores whether the last getNextVariant call succeeded. If it didn't
    // succeed, we can't call it again (without setRegion-ing) or vcflib will
    // crash.
    bool safe_to_get = true;
    // This is the actual buffer.
    vcflib::Variant buffer;
    // This points to the VariantCallFile we wrap
    // We can wrap the null file (and never have any variants) with a null here.
    vcflib::VariantCallFile* const file;
    
    

private:
    // Don't copy these or it will break the buffering semantics.
    VcfBuffer(const VcfBuffer& other) = delete;
    VcfBuffer& operator=(const VcfBuffer& other) = delete;
        
};

}

#endif
