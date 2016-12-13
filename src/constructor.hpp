#ifndef VG_CONSTRUCTOR_HPP
#define VG_CONSTRUCTOR_HPP

/** \file
 * constructor.hpp: defines a tool class used for constructing VG graphs from
 * VCF files.
 */

#include <vector>
#include <set>
#include <map>
#include <cstdlib>
#include <functional>
#include <regex>

#include "types.hpp"
#include "progressive.hpp"

#include "vg.pb.h"

// We need vcflib
#include "Variant.h"
// And fastahack
#include "Fasta.h"

namespace vg {

using namespace std;

/**
 * Represents a constructed region of the graph alogn a single linear sequence.
 * Contains the protobuf Graph holding all the created components (which may be
 * too large to serialize), a set of node IDs whose left sides need to be
 * connected to when you connect to the start of the chunk, and a set of node
 * IDs whose right sides need to be connected to when you connect to the end of
 * the chunk.
 *
 * Node ordering is restricted: if there is a single source, it must be the very
 * first node in the graph with ID 1, and if there is a single sink it must be
 * the very last node in the graph with ID max_id. Additionally, single sources
 * and single sinks must be visited by only a single path, the reference path.
 *
 * The overall reference path must also always be path 0. Also, all mappings in
 * all paths must be full-length matches on the forward strand, and they must be
 * sorted by rank. Ranks must be filled and start with rank 1 in each path.
 */
struct ConstructedChunk {
    // What nodes, edges, and mappings exist?
    // The graph will start at node 1.
    Graph graph;
    
    // What's the max ID used? This is useful for bumping other chunks out of
    // the way.
    id_t max_id;
    
    // What nodes have left sides that match up with the left edge of the chunk?
    set<id_t> left_ends;
    // And similarly for right sides on the right edge of the chunk?
    set<id_t> right_ends;
};

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

class Constructor : public Progressive {

public:

    // Should alts be interpreted as flat (false) or aligned back to the
    // reference by vcflib (true)?
    bool flat = false;
    
    // Should we add paths for the different alts of variants, like
    // _alt_6079b4a76d0ddd6b4b44aeb14d738509e266961c_0 and
    // _alt_6079b4a76d0ddd6b4b44aeb14d738509e266961c_1?
    bool alt_paths = false;

    // Should we handle structural variants in the VCF file,
    // or at least the ones we know how to?
    bool do_svs = false;

    // Should we also store the alt_paths as loci?
    // e.g.
    // Locus{
    //  Name: locus1,
    //  paths: [alt1, alt2, alt3]
    // }
    //
    // Hacky variant index: map<string, Locus> where
    // string is a string name for the VCF entry (referenceID_contig_pos_SVTYPE_hash(alt_sequence))
    bool alts_as_loci = false;
    
    // If true, break boring sequence into pieces greedily. If false, divide
    // over-long nodes into more even pieces.
    bool greedy_pieces = false;

    // If a deletion deletes the anchoring base of another deletion (or just
    // occurs directly adjacent to another deletion), can we chain them in the
    // graph to allow the longer combined deletion?
    bool chain_deletions = true;
    
    // What's the maximum node size we should allow?
    size_t max_node_size = 1000;
    
    // How many variants do we want to put into a chunk? We'll still go over
    // this by a bit when we fetch all the overlapping variants, but this is how
    // many we shoot for.
    size_t vars_per_chunk = 1024;
    
    // How many bases do we want to have per chunk? We don't necessarily want to
    // load all of chr1 into an std::string, even if we have no variants on it.
    size_t bases_per_chunk = 1024 * 1024;
    
    /**
     * Add a name mapping between a VCF contig name and a FASTA sequence name.
     * Both must be unique.
     */
    void add_name_mapping(const string& vcf_name, const string& fasta_name);
    
    /**
     * Convert the given VCF contig name to a FASTA sequence name, through the
     * rename mappings.
     */
    string vcf_to_fasta(const string& vcf_name) const;
    
    /**
     * Convert the given FASTA sequence name to a VCF contig name, through the
     * rename mappings.
     */
    string fasta_to_vcf(const string& fasta_name) const;
    
    // This set contains the set of VCF sequence names we want to build the
    // graph for. If empty, we will build the graph for all sequences in the
    // FASTA. If nonempty, we build only for the specified sequences. If
    // vcf_renames applies a translation, these should be pre-translation, VCF-
    // namespace names.
    set<string> allowed_vcf_names;
    
    // This map maps from VCF sequence name to a (start, end) interval, in
    // 0-based end-exclusive coordinates, for the region of the sequence to
    // include in the graph. If it is set for a sequence, only that part of the
    // VCF will be used, and only that part of the primary path will be present
    // in the graph. If it is unset, the whole contig's graph will be
    // constructed. If vcf_renames applies a translation, keys should be pre-
    // translation, VCF-namespace names.
    map<string, pair<size_t, size_t>> allowed_vcf_regions;
    
    /**
     * Construct a ConstructedChunk of graph from the given piece of sequence,
     * with the given name, applying the given variants. The variants need to be
     * sorted by start position, and have their start positions set to be ZERO-
     * BASED. However, they also need to have their start positions relative to
     * the global start of the contig, so that hash-based names come out right
     * for them. They also need to not overlap with any variants not in the
     * vector we have (i.e. we need access to all overlapping variants for this
     * region). The variants must not extend beyond the given sequence, though
     * they can abut its edges.
     *
     * Variants in the vector may not use symbolic alleles.
     *
     * chunk_offset gives the global 0-based position at which this chunk starts
     * in the reference contig it is part of, which is used to correctly place
     * variants.
     */
    ConstructedChunk construct_chunk(string reference_sequence, string reference_path_name,
        vector<vcflib::Variant> variants, size_t chunk_offset) const;
    
    /**
     * Construct a graph for the given VCF contig name, using the given
     * reference and the variants from the given buffered VCF file. Emits a
     * sequence of Graph chunks, which may be too big to serealize directly.
     *
     * Doesn't handle any of the setup for VCF indexing. Just scans all the
     * variants that can come out of the buffer, so make sure indexing is set on
     * the file first before passing it in.
     */
    void construct_graph(string vcf_contig, FastaReference& reference, VcfBuffer& variant_source,
         const vector<FastaReference*>& insertion, function<void(Graph&)> callback);
    
    /**
     * Construct a graph using the given FASTA references and VCFlib VCF files.
     * The VCF files are assumed to be grouped by contig and then sorted by
     * position within the contig, such that each contig is present in only one
     * file. If multiple FASTAs are used, each contig must be present in only
     * one FASTA file. Reference and VCF vectors may not contain nulls.
     */
    void construct_graph(const vector<FastaReference*>& references, const vector<vcflib::VariantCallFile*>& variant_files,
        const vector<FastaReference*>& insertions, function<void(Graph&)> callback);
    
protected:
    
    /// This map maps from VCF sequence names to FASTA sequence names. If a
    /// VCF sequence name doesn't appear in here, it gets passed through
    /// unchanged. Note that the primary path for each contig will be named after
    /// the FASTA sequence name and not the VCF sequence name.
    map<string, string> vcf_to_fasta_renames;
    
    /// This is the reverse map from FASTA sequence name to VCF sequence name.
    map<string, string> fasta_to_vcf_renames;

    /// Remembers which unusable symbolic alleles we've already emitted a warning
    /// about during construction.
    set<string> symbolic_allele_warnings;

private:

    /**
     * Given a vector of lists of VariantAllele edits, trim in from the left and
     * right, leaving a core of edits bounded by edits that actually change the
     * reference in at least one allele.
     *
     * Postcondition: either all lists of VariantAlleles are empty, or at least
     * one begins with a non-match and at least one ends with a non-match.
     */
    static void trim_to_variable(vector<list<vcflib::VariantAllele>>& parsed_alleles);
    
    /**
     * Given a list of VariantAllele edits, condense adjacent perfect match
     * edits to be maximally long.
     */
    static void condense_edits(list<vcflib::VariantAllele>& parsed_allele); 

    /**
     * Given a vector of lists of VariantAllele edits that have been trimmed
     * with trim_to_variable() above, one per non-reference alt for a variant,
     * return the position of the first varaible base, and the position of the
     * last variable base. If there's no variable-region, the result is max
     * int64_t and -1, and if there's a 0-length variable region, the result is
     * the base after it and the base before it.
     */
    static pair<int64_t, int64_t> get_bounds(const vector<list<vcflib::VariantAllele>>& trimmed_variant);
    

};

}

#endif
