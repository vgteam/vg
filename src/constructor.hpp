#ifndef VG_CONSTRUCTOR_HPP_INCLUDED
#define VG_CONSTRUCTOR_HPP_INCLUDED

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
#include "vcf_buffer.hpp"
#include "name_mapper.hpp"
#include "handle.hpp"

#include <vg/vg.pb.h>

// We need vcflib
#include <Variant.h>
// And fastahack
#include <Fasta.h>



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

class Constructor : public Progressive, public NameMapper {

public:

    // Should alts be interpreted as flat (false) or aligned back to the
    // reference by vcflib (true)?
    bool flat = false;
    
    // In non-flat mode, how big can the longest allele of a variant be before
    // we fall back to flat mode anyway?
    size_t max_parsed_variant_size = 100;
    
    // Should we add paths for the different alts of variants, like
    // _alt_6079b4a76d0ddd6b4b44aeb14d738509e266961c_0 and
    // _alt_6079b4a76d0ddd6b4b44aeb14d738509e266961c_1?
    bool alt_paths = false;

    // Should we handle structural variants in the VCF file,
    // or at least the ones we know how to?
    bool do_svs = false;

    // Should we trim the 1bp reference sequence that by default is placed
    // on indel variants?
    bool trim_indels = true;

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
    
    // Should we warn if lowercase characters are encountered in each input sequence?
    bool warn_on_lowercase = true;
    
    // Should we warn if IUPAC ambiguity codes (other than N) are encountered
    // in each input sequence?
    bool warn_on_ambiguous = true;
    
    // What's the maximum node size we should allow?
    size_t max_node_size = 1000;
    
    // How many variants do we want to put into a chunk? We'll still go over
    // this by a bit when we fetch all the overlapping variants, but this is how
    // many we shoot for.
    size_t vars_per_chunk = 1024;
    
    // How many bases do we want to have per chunk? We don't necessarily want to
    // load all of chr1 into an std::string, even if we have no variants on it.
    size_t bases_per_chunk = 1024 * 1024;
    
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
     *
     * insertion contains FASTAs containing serquences for resolving symbolic
     * insert alleles in the VCF.
     *
     * Calls the given callback with constructed graph chunks, in a single
     * thread. Chunks may contain dangling edges into the next chunk.
     */
    void construct_graph(string vcf_contig, FastaReference& reference, VcfBuffer& variant_source,
         const vector<FastaReference*>& insertion, const function<void(Graph&)>& callback);
    
    /**
     * Construct a graph using the given FASTA references and VCFlib VCF files.
     * The VCF files are assumed to be grouped by contig and then sorted by
     * position within the contig, such that each contig is present in only one
     * file. If multiple FASTAs are used, each contig must be present in only
     * one FASTA file. Reference and VCF vectors may not contain nulls.
     *
     * insertions contains FASTAs containing serquences for resolving symbolic
     * insert alleles in the VCFs.
     *
     * Calls the given callback with constructed graph chunks, eventually
     * (hopefully) in multiple threads. Chunks may contain dangling edges into
     * the next chunk.
     */
    void construct_graph(const vector<FastaReference*>& references, const vector<vcflib::VariantCallFile*>& variant_files,
        const vector<FastaReference*>& insertions, const function<void(Graph&)>& callback);
        
    /**
     * Construct a graph using the given FASTA references and VCF files on disk.
     * The VCF files are assumed to be grouped by contig and then sorted by
     * position within the contig, such that each contig is present in only one
     * file. If multiple FASTAs are used, each contig must be present in only
     * one FASTA file.
     *
     * insertions contains FASTA filenames containing serquences for resolving
     * symbolic insert alleles in the VCFs.
     *
     * Calls the given callback with constructed graph chunks, eventually
     * (hopefully) in multiple threads. Chunks may contain dangling edges into
     * the next chunk.
     */
    void construct_graph(const vector<string>& reference_filenames, const vector<string>& variant_filenames,
        const vector<string>& insertion_filenames, const function<void(Graph&)>& callback);
        
    /**
     * Construct a graph using the given FASTA references and VCFlib VCF files.
     * The VCF files are assumed to be grouped by contig and then sorted by
     * position within the contig, such that each contig is present in only one
     * file. If multiple FASTAs are used, each contig must be present in only
     * one FASTA file. Reference and VCF vectors may not contain nulls.
     *
     * insertions contains FASTAs containing serquences for resolving symbolic
     * insert alleles in the VCFs.
     *
     * Builds the graph into the given mutable graph object, which may not be
     * thread safe.
     */
    void construct_graph(const vector<FastaReference*>& references, const vector<vcflib::VariantCallFile*>& variant_files,
        const vector<FastaReference*>& insertions, MutablePathMutableHandleGraph* destination);
        
    /**
     * Construct a graph using the given FASTA references and VCF files on disk.
     * The VCF files are assumed to be grouped by contig and then sorted by
     * position within the contig, such that each contig is present in only one
     * file. If multiple FASTAs are used, each contig must be present in only
     * one FASTA file.
     *
     * insertions contains FASTA filenames containing serquences for resolving
     * symbolic insert alleles in the VCFs.
     *
     * Builds the graph into the given mutable graph object, which may not be
     * thread safe.
     */
    void construct_graph(const vector<string>& reference_filenames, const vector<string>& variant_filenames,
        const vector<string>& insertion_filenames, MutablePathMutableHandleGraph* destination);
    
protected:
    
    /// Remembers which unusable symbolic alleles we've already emitted a warning
    /// about during construction.
    set<string> symbolic_allele_warnings;
    
    /// All chunks are generated with IDs starting at 1, but graphs emitted from
    /// construct_graph need to have the IDs rewritten so they don't overlap.
    /// Moreover, multiple calls to construct_graph need to not have conflicting
    /// IDs, because some construct_graph implementations call other ones. What
    /// we do for now is globally track the max ID already used, so all calls to
    /// construct_graph follow a single ID ordering.
    id_t max_id = 0;

private:

    /**
     * Given a vector of lists of VariantAllele edits, trim in from the left and
     * right, leaving a core of edits bounded by edits that actually change the
     * reference in at least one allele.
     *
     * Postcondition: either all lists of VariantAlleles are empty, or at least
     * one begins with a non-match and at least one ends with a non-match.
     * Adjacent edits in the list abut; there are no uncovered gaps in the edits.
     * This means that *internal* perfect match edits will be preserved.
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
    
    /**
     * Given a symbolic variant, check its bounds and return them. This
     * function is needed to handle SVs properly, since they won't always have
     * their ref and alt fields put in. Note that insertions may have an end
     * bound before their start, because the anchoring base isn't included.
     */
    static pair<int64_t, int64_t> get_symbolic_bounds(vcflib::Variant var);
    
    /**
     * Given a sequence, get rid of all the lowercase characters and all the
     * ambiguity codes. Warn if configured, and the sequence has a name
     * assigned, and no warning has yet been issued for that name, or if a
     * variant is specified.
     *
     * Will error if this results in a string with anything other than A, C, G,
     * T, and N.
     *
     * sequence_start_offset can be set to produce useful messages if the
     * sequence we are looking at is an excerpt from a longer sequence.
     *
     * Santitizing may move the stored string data in memory.
     *
     * Returns true if the string was modified.
     *
     * We need this as a function because vcflib reaches back and reads the
     * FASTA files directly, so we can't *just* preprocess the reference and we
     * need to constantly clean up the variants.
     */
    bool sanitize_sequence_in_place(string& sequence, const string* sequence_name = nullptr, size_t sequence_start_offset = 0, const vcflib::Variant* variant = nullptr) const;
    
    /// What sequences have we warned about containing lowercase characters?
    mutable unordered_set<string> lowercase_warned_sequences;
    /// Have we given a warning yet about lowercase alt alleles?
    mutable bool lowercase_warned_alt = false;
    /// Have we given a warning yet about multiallelic SVs?
    mutable bool multiallelic_sv_warned = false;
    /// Have we given a warning yet about uncanonicalizable SVs?
    mutable bool uncanonicalizable_sv_warned = false;
    /// What sequences have we warned about containing unsupported ambiguity codes?
    mutable unordered_set<string> ambiguous_warned_sequences;
    

};

}

#endif
