#ifndef VG_VCF_BUFFER_HPP_INCLUDED
#define VG_VCF_BUFFER_HPP_INCLUDED

/** \file
 * vcf_buffer.hpp: defines a buffer for reading through VCF files with
 * look-ahead.
 */

#include <list>
#include <tuple>

// We need vcflib
#include "Variant.h"


namespace vg {

using namespace std;

/**
 * Provides a one-variant look-ahead buffer on a vcflib::VariantFile. Lets
 * construction functions peek and see if they want the next variant, or lets
 * them ignore it for the next construction function for a different contig to
 * handle. Ought not to be copied.
 */
class VcfBuffer {

public:
    /**
     * Return a pointer to the buffered variant, or null if no variant is
     * buffered. Pointer is invalidated when the buffer is handled.
     *
     * Although the variant is not const, it may not be moved out of or modified
     * in ways that confuse vcflib.
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
    bool has_tabix() const;
    
    /**
     * This tries to set the region on the underlying vcflib VariantCallFile to
     * the given contig and region, if specified. Coordinates coming in should
     * be 0-based, and will be converted to 1-based internally.
     *
     * Returns true if the region was successfully set, and false otherwise (for
     * example, if there is not tabix index, or if the given region is not part
     * of this VCF. Note that if there is a tabix index, and set_region returns
     * false, the position in the VCF file is undefined until the next
     * successful set_region call.
     *
     * If either of start and end are specified, then both of start and end must
     * be specified.
     *
     * Discards any variants previously in the buffer.
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

/**
 * Provides a look-around buffer for VCFs where you can look at each variant in
 * the context of nearby variants.
 *
 * Also caches parsings of genotypes, so you can iterate over genotypes
 * efficiently without parsing them out over and over again.
 */
class WindowedVcfBuffer {

public:

    /**
     * Make a new WindowedVcfBuffer buffering the file at the given pointer
     * (which must outlive the buffer, but which may be null). The VCF in the
     * file must be sorted, but may contain overlapping variants.
     */
    WindowedVcfBuffer(vcflib::VariantCallFile* file, size_t window_size);
    
    /**
     * Advance to the next variant, making it the current variant. Returns true
     * if a next variant exists, and false if no next variant can be found. Must
     * be called (and return true) before the first call to get() after
     * constructing the WindowedVcfBuffer or setting the region.
     */
    bool next();
    
    /**
     * Get the current variant in its context. Throws an exception if no variant
     * is current. Returns a vector of variants in the window before the current
     * variant, the current variant, and a vector of variants in the window
     * after the current variant.
     *
     * Pointers will be invalidated upon the next call to next() or
     * set_region().
     */
    tuple<vector<vcflib::Variant*>, vcflib::Variant*, vector<vcflib::Variant*>> get();
    
    /**
     * Like get(), but elides variants in the context that overlap the current
     * variant, or each other.
     */
    tuple<vector<vcflib::Variant*>, vcflib::Variant*, vector<vcflib::Variant*>> get_nonoverlapping();
    
    /**
     * Given a pointer to a cached variant owned by this WindowedVcfBuffer (such
     * as might be obtained from get()), return the cached parsed-out genotypes
     * for all the samples, in the order the samples appear in the VCF file.
     *
     * Returns a reference which is valid until the variant passed in is
     * scrolled out of the buffer.
     */
    const vector<vector<int>>& get_parsed_genotypes(vcflib::Variant* variant);
    
    /**
     * This returns true if we have a tabix index, and false otherwise. If this
     * is false, set_region may be called, but will do nothing and return false.
     */
    bool has_tabix() const;
    
    /**
     * This tries to set the region on the underlying vcflib VariantCallFile to
     * the given contig and region, if specified. Coordinates coming in should
     * be 0-based, and will be converted to 1-based internally.
     *
     * Returns true if the region was successfully set, and false otherwise (for
     * example, if there is not tabix index, or if the given region is not part
     * of this VCF. Note that if there is a tabix index, and set_region returns
     * false, the position in the VCF file is undefined until the next
     * successful set_region call.
     *
     * If either of start and end are specified, then both of start and end must
     * be specified.
     *
     * Discards any variants previously in the buffer.
     */
    bool set_region(const string& contig, int64_t start = -1, int64_t end = -1);
    
protected:
    
    /**
     * Quickly decompose a genotype without any string copies.
     */
    static vector<int> decompose_genotype_fast(const string& genotype);
    
    // This lets us read from our VCF
    VcfBuffer reader;
    
    // This holds the window size around the start of the current variant, in
    // bp. We are only interested in variants that start in this window.
    size_t window_size;
    
    // This holds all the variants within the window before the current variant
    list<unique_ptr<vcflib::Variant>> variants_before;
    // This holds all the variants within the window after the current variant
    list<unique_ptr<vcflib::Variant>> variants_after;
    
    // The variant we are currently "on", which we are notionally in the process
    // of moving from the after list to the before list.
    unique_ptr<vcflib::Variant> current;
    
    // Keep a cache of parsed genotypes for variants that are currently in one
    // of our lists. Genotypes are vectors of allele numbers, one per sample,
    // and occur in the order that the samples occur in the file.
    map<vcflib::Variant*, vector<vector<int>>> cached_genotypes;
    
    // Keep a key from sample name index in map key order to sample index in the
    // VCF. We use this to build our cache efficiently without any lample name
    // lookups.
    vector<size_t> map_order_to_original;
    
private:
    // Don't copy or assign because we contain VcfBuffers
    WindowedVcfBuffer(const WindowedVcfBuffer& other) = delete;
    WindowedVcfBuffer& operator=(const WindowedVcfBuffer& other) = delete;

};

}

#endif
