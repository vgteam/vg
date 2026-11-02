#include <atomic>
#include <charconv>
#include <chrono>
#include <cstdio>
#include <limits>

#include <omp.h>

#include "sorted_vcf_writer.hpp"
#include "graph_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"
#include "utility.hpp"

//#define debug

namespace vg {

bool SortedVCFWriter::BufferedRecordKey::operator<(const BufferedRecordKey& other) const {
    return std::tie(contig, position, id, block)
           < std::tie(other.contig, other.position, other.id, other.block);
}

void SortedVCFWriter::add_info_line(string line) {
    info_lines.push_back(std::move(line));
}

void SortedVCFWriter::add_format_line(string line) {
    format_lines.push_back(std::move(line));
}

void SortedVCFWriter::add_filter_line(string line) {
    filter_lines.push_back(std::move(line));
}

const NodeTranslation* SortedVCFWriter::translation() const {
    return node_translation;
}

const string& SortedVCFWriter::sample_name() const {
    return sample;
}

Region record_ref_region(const vcflib::Variant& record) {
    return Region{record.sequenceName, record.position,
                  record.position + (int64_t)record.ref.length() - 1};
}

SortedVCFWriter::SortedVCFWriter(const string& sample_name) : sample(sample_name) {
    output_vcf.sampleNames = {sample_name};
    output_variants.resize(get_thread_count());
    suppressed_ref_info.resize(get_thread_count());
}

void SortedVCFWriter::set_sample_names(vector<string> names) {
    output_vcf.sampleNames = std::move(names);
}

void SortedVCFWriter::set_contigs(vector<ReferenceContig> declared, bool only_with_records) {
    contigs = std::move(declared);
    only_contigs_with_records = only_with_records;
}

string SortedVCFWriter::vcf_header(const unordered_set<string>* keep_contigs) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;
    for (const ReferenceContig& contig : contigs) {
        if (keep_contigs != nullptr && !keep_contigs->contains(contig.name)) {
            // A reference contig that has no record is not worth declaring: with a gref cover
            // most contigs are fragments, and on a human chromosome a third of them carry
            // nothing.
            continue;
        }
        ss << "##contig=<ID=" << contig.name << ",length=" << contig.length << ">" << endl;
    }
    if (include_nested) {
        ss << nesting_info_headers();
    }
    for (const string& line : info_lines) {
        ss << line << endl;
    }
    ss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    for (const string& line : format_lines) {
        ss << line << endl;
    }
    ss << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
    for (const string& line : filter_lines) {
        ss << line << endl;
    }
    for (const string& name : output_vcf.sampleNames) {
        ss << "##SAMPLE=<ID=" << name << ">" << endl;
    }
    ss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const string& name : output_vcf.sampleNames) {
        ss << "\t" << name;
    }
    ss << endl;
    return ss.str();
}


bool SortedVCFWriter::add_variant(vcflib::Variant& var, size_t block, bool hidden) {
    if (hidden) {
        // Only the nesting tags read a hidden record.
        if (include_nested) {
            suppressed_ref_info[omp_get_thread_num()][var.id] = record_ref_region(var);
        }
        return true;
    }
    var.setVariantCallFile(output_vcf);
    stringstream ss;
    ss << var;
    string dest;
    if (ss.str().length() > SortedVCFWriter::max_vcf_line_length) {
        return false;
    }         
    int ret = zstdutil::CompressString(ss.str(), dest);
    assert(ret == 0);
    // the Variant object is too big to keep in memory when there are many genotypes, so we
    // store it in a zstd-compressed string
    output_variants[omp_get_thread_num()].emplace_back(
        BufferedRecordKey{var.sequenceName, (size_t)var.position, var.id, block}, dest);
    return true;
}

void SortedVCFWriter::write_variants(
    ostream& out_stream, const SnarlDecomposition* sites, const HandleGraph* graph,
    const function<string(const net_handle_t&)>& record_id_of_site) {
    // Find the contigs with records now, since sorting empties the buffers.
    unordered_set<string> output_contigs;
    if (only_contigs_with_records) {
        output_contigs = get_output_contigs();
    }
    string header = vcf_header(only_contigs_with_records ? &output_contigs : nullptr);
    // The file gets a copy of the header, because parsing removes the last newline.
    string parsed_header = header;
    output_vcf.openForOutput(parsed_header);
    out_stream << header;
    if (include_nested) {
        assert(sites != nullptr);
        if (record_id_of_site) {
            update_nesting_info_tags(*sites, record_id_of_site);
        } else {
            // Name each snarl from its boundary nodes, read the way the snarl is traversed.
            assert(graph != nullptr);
            update_nesting_info_tags(*sites, [&](const net_handle_t& net) {
                bool backward = sites->starts_at(net) == SnarlDecomposition::END;
                handle_t start = sites->get_handle(sites->get_bound(net, backward, true), graph);
                handle_t end = sites->get_handle(sites->get_bound(net, !backward, false), graph);
                return site_name(*graph, start, end, node_translation, false);
            });
        }
    }
    vector<pair<BufferedRecordKey, string>> all_variants;
    // Reserve room for every record at once, rather than growing for each buffer.
    size_t total_variants = 0;
    for (const auto& buf : output_variants) {
        total_variants += buf.size();
    }
    all_variants.reserve(total_variants);
    // Move the records out of each buffer and free it. `buf` is not const, because moving from a
    // const range copies, and the whole VCF is in memory here.
    for (auto& buf : output_variants) {
        std::ranges::move(buf, std::back_inserter(all_variants));
        buf.clear();
        buf.shrink_to_fit();
    }
    std::sort(all_variants.begin(), all_variants.end(),
              [](const pair<BufferedRecordKey, string>& v1,
                 const pair<BufferedRecordKey, string>& v2) {
                  return v1.first < v2.first;
              });

    // Decompress the records in parallel, a batch at a time, and write each batch in order. Only
    // one batch of text is in memory at once.
    const size_t batch_records = 1 << 16;
    vector<string> lines;
    for (size_t batch_start = 0; batch_start < all_variants.size(); batch_start += batch_records) {
        const size_t batch_end = min(all_variants.size(), batch_start + batch_records);
        lines.assign(batch_end - batch_start, string());
#pragma omp parallel for schedule(dynamic, 256)
        for (size_t record_i = batch_start; record_i < batch_end; ++record_i) {
            const string& compressed = all_variants[record_i].second;
            string& dest = lines[record_i - batch_start];
            int ret = zstdutil::DecompressString(compressed, dest);
            assert(ret == 0);
        }
        for (const string& line : lines) {
            // endl would flush, making one write per record.
            out_stream << line << '\n';
        }
    }
}



static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                if (g->first == alleleIndex) {
                    alts += g->second;
                }
            }
        }
    }
    return alts;
}

static int countAlleles(vcflib::Variant& var) {
    int alleles = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
		    alleles += g->second;
		}
            }
        }
    }
    return alleles;
}

// this isn't from vcflib, but seems to make more sense than just returning the number of samples in the file again and again
static int countSamplesWithData(vcflib::Variant& var) {
    int samples_with_data = 0;
    for (auto s = var.samples.begin(); s != var.samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator gt = sample.find("GT");
        bool has_data = false;
        if (gt != sample.end()) {
            map<int, int> genotype = vcflib::decomposeGenotype(gt->second.front());
            for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
		if (g->first != vcflib::NULL_ALLELE) {
                    has_data = true;
                    break;
		}
            }
        }
        if (has_data) {
            ++samples_with_data;
        }
    }
    return samples_with_data;
}

void SortedVCFWriter::vcf_fixup(vcflib::Variant& var) const {
    // copied from https://github.com/vgteam/vcflib/blob/master/src/vcffixup.cpp
    
    stringstream ns;
    ns << countSamplesWithData(var);
    var.info["NS"].clear();
    var.info["NS"].push_back(ns.str());

    var.info["AC"].clear();
    var.info["AF"].clear();
    var.info["AN"].clear();

    int allelecount = countAlleles(var);
    stringstream an;
    an << allelecount;
    var.info["AN"].push_back(an.str());

    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string& allele = *a;
        int altcount = countAlts(var, var.getAltAlleleIndex(allele) + 1);
        stringstream ac;
        ac << altcount;
        var.info["AC"].push_back(ac.str());
        stringstream af;
        double faf = (double) altcount / (double) allelecount;
        if(faf != faf) faf = 0;
        af << faf;
        var.info["AF"].push_back(af.str());
    }
}

void SortedVCFWriter::set_translation(const NodeTranslation* translation) {
    node_translation = translation;
}

void SortedVCFWriter::set_nested(bool nested) {
    include_nested = nested;
}

void SortedVCFWriter::set_gref_levels(map<string, int> levels) {
    this->gref_levels = std::move(levels);
}

unordered_set<string> SortedVCFWriter::get_output_contigs() const {
    unordered_set<string> with_records;
    // The contig is in each record's sort key, so nothing has to be decompressed.
    for (const auto& thread_buf : output_variants) {
        for (const auto& [key, record] : thread_buf) {
            with_records.insert(key.contig);
        }
    }
    return with_records;
}






string SortedVCFWriter::nesting_info_headers() {
    stringstream ss;
    ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree counting only ancestors whose record is on this record's own reference contig (0=top level for this contig)\">" << endl;
    ss << "##INFO=<ID=CH,Number=1,Type=Integer,Description=\"Nesting steps between VCF reference contigs: how many coordinate-system changes separate this record from a linear reference. Counted as the greater of the in-VCF ancestor hops and the record's own gref contig level, because counting only ancestors that happened to emit a record made a record on a gref fragment whose parent produced no line indistinguishable from one on the linear reference -- 29,843 of 41,669 off-reference records on a gref-covered chr20. So CH >= 1 no longer implies an in-VCF parent, and therefore no longer implies PS\">" << endl;
    ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    ss << "##INFO=<ID=RC,Number=1,Type=String,Description=\"CHROM of the topmost ancestor record in this VCF, or this record's own CHROM when it has none. On a gref fragment, where that own CHROM would be no use, the enclosing site is named even if it produced no record of its own; the tags are absent when there is no such site either\">" << endl;
    ss << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Start of the site named by RC: the POS of its record, or where the site begins when it produced none. A position on that contig, not a span of the snarl, so it can precede the sequence this record describes\">" << endl;
    ss << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"End of the site named by RC: RS plus the length of that site's REF allele\">" << endl;
    return ss.str();
}


namespace {

/// Return true if `net` is a site: a snarl, or an element that is none of the kinds a
/// SnarlDecomposition knows, such as an edit in an EditSnarlDecompositionOverlay.
bool is_site(const SnarlDecomposition& sites, const net_handle_t& net) {
    return sites.is_snarl(net)
           || !(sites.is_chain(net) || sites.is_node(net) || sites.is_sentinel(net)
                || sites.is_root(net));
}

/// Return the site that `net` lies directly inside, or a root handle if none does.
net_handle_t parent_site(const SnarlDecomposition& sites, const net_handle_t& net) {
    net_handle_t here = sites.get_parent(net);
    while (!sites.is_root(here) && !sites.is_snarl(here)) {
        here = sites.get_parent(here);
    }
    return here;
}

/// Call `visit` on `top`, if it is a site, and on each site inside it, each before the sites
/// inside it.
///
/// A site other than a snarl has no children. Any other element, even a node, can have a site as
/// a child, so the search goes into every one.
void for_each_site_in(const SnarlDecomposition& sites, const net_handle_t& top,
                      const function<void(const net_handle_t&)>& visit) {
    // The elements left to search, the next one last.
    vector<net_handle_t> stack{top};
    vector<net_handle_t> children;
    while (!stack.empty()) {
        net_handle_t here = stack.back();
        stack.pop_back();
        if (is_site(sites, here)) {
            visit(here);
            if (!sites.is_snarl(here)) {
                continue;
            }
        }
        children.clear();
        sites.for_each_child(here, [&](const net_handle_t& child) {
            children.push_back(child);
        });
        stack.insert(stack.end(), children.rbegin(), children.rend());
    }
}

/// Call `visit` on every site of `sites`.
///
/// With `in_preorder`, call it on this thread, on each site before the sites inside it. Otherwise
/// call it on several threads, one child of the root at a time, in any order.
void for_each_site(const SnarlDecomposition& sites,
                   const function<void(const net_handle_t&)>& visit, bool in_preorder) {
    vector<net_handle_t> tops;
    sites.for_each_child(sites.get_root(), [&](const net_handle_t& child) {
        tops.push_back(child);
    });
    if (in_preorder) {
        for (const net_handle_t& top : tops) {
            for_each_site_in(sites, top, visit);
        }
        return;
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < tops.size(); ++i) {
        for_each_site_in(sites, tops[i], visit);
    }
}

/// Map each record name in `chrom_of_name` to the site that has it.
///
/// A site's record can be written in either orientation, so each site is looked up by the ID of
/// its record read both ways round. Where two sites have one name, the name goes to the one that
/// the preorder visits last.
unordered_map<string, net_handle_t>
index_sites_by_record_name(const SnarlDecomposition& sites,
                           const unordered_map<string, uint32_t>& chrom_of_name,
                           const function<string(const net_handle_t&)>& record_id_of_site) {
    auto for_each_name = [&](const net_handle_t& site, const auto& match) {
        string name = record_id_of_site(site);
        string flipped_name = record_id_of_site(sites.flip(site));
        if (auto named = chrom_of_name.find(name); named != chrom_of_name.end()) {
            match(named->first);
        }
        if (flipped_name != name) {
            if (auto named = chrom_of_name.find(flipped_name); named != chrom_of_name.end()) {
                match(named->first);
            }
        }
    };
    // Match the sites on several threads, each into a list of its own. The matches depend on the
    // order of the sites only if two different sites match one name, and then we match again in
    // preorder.
    vector<vector<pair<const string*, net_handle_t>>> found(max(1, omp_get_max_threads()));
    for_each_site(sites, [&](const net_handle_t& site) {
        auto& mine = found[omp_get_thread_num()];
        for_each_name(site, [&](const string& name) { mine.emplace_back(&name, site); });
    }, false);
    unordered_map<string, net_handle_t> name_to_snarl;
    name_to_snarl.reserve(chrom_of_name.size());
    bool ambiguous = false;
    for (const auto& thread_found : found) {
        for (const auto& [name, site] : thread_found) {
            auto [placed, inserted] = name_to_snarl.try_emplace(*name, site);
            if (!inserted && placed->second != site) {
                ambiguous = true;
            }
        }
    }
    if (ambiguous) {
        name_to_snarl.clear();
        for_each_site(sites, [&](const net_handle_t& site) {
            for_each_name(site, [&](const string& name) { name_to_snarl[name] = site; });
        }, true);
    }
    return name_to_snarl;
}

/// Read the name of the site of each buffered record, and the stretch of the reference its REF
/// covers, decompressing the buffers in parallel. The name is the same for every record of a site
/// split into several.
///
/// The result has one list per buffer, in the order of the buffer.
vector<vector<pair<string, Region>>> read_record_regions(
    const vector<vector<pair<SortedVCFWriter::BufferedRecordKey, string>>>& buffers) {
    vector<vector<pair<string, Region>>> placed(buffers.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < buffers.size(); ++b) {
        placed[b].reserve(buffers[b].size());
        string output_variant_string;
        for (const auto& [key, compressed] : buffers[b]) {
            output_variant_string.clear();
            int ret = zstdutil::DecompressString(compressed, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            int64_t position = stoll(toks[1]);
            placed[b].emplace_back(block_site_name(toks[2]),
                                   Region{std::move(toks[0]), position,
                                          position + (int64_t)toks[3].length() - 1});
        }
    }
    return placed;
}

}

void SortedVCFWriter::update_nesting_info_tags(
    const SnarlDecomposition& sites,
    const function<string(const net_handle_t&)>& record_id_of_site) {

    auto name_of = [&](const net_handle_t& site) { return record_id_of_site(site); };
    auto flipped_name_of = [&](const net_handle_t& site) {
        return record_id_of_site(sites.flip(site));
    };

    // Merge the per-thread suppressed-site intervals collected during calling.  These are sites
    // that never reached the VCF, so pass 1 below cannot see them, but a record nested under one
    // has no other way to name a reference position.
    unordered_map<string, Region> suppressed_ref;
    for (auto& buf : suppressed_ref_info) {
        for (auto& kv : buf) {
            suppressed_ref.emplace(kv.first, std::move(kv.second));
        }
        buf.clear();
        buf.rehash(0);
    }

    // pass 1) index sites in vcf
    // (todo: this could be done more quickly upstream)
    //
    // One index, not two: presence in chrom_of_name IS "this snarl name is in the VCF", and
    // the value is which reference contig its record landed on.  Keeping a separate
    // names_in_vcf set alongside would store all 400k snarl-ID strings twice, which measured
    // as +70 MB of peak RSS on chr22 -- the keys, not the values, are what costs.
    //
    // Contig names are interned rather than stored per record for the same reason: there are
    // at most a few thousand distinct ones, and all we ever ask is whether two are the same.
    unordered_map<string, uint32_t> chrom_index;
    // Whether each interned contig is a synthetic gref fragment, by the same index.  Stored as a
    // bit per contig rather than looked up by name later, so the names are still stored once.
    // is_gref_name(), not is_gref_derived(): a gref copy of a real reference contig is a perfectly
    // good coordinate system -- it is what the whole VCF is deconstructed against -- and only the
    // "_<N>_alt" fragments are positions a reader cannot look up.
    vector<bool> chrom_is_gref_fragment;
    auto intern_chrom = [&](const string& chrom) -> uint32_t {
        auto result = chrom_index.emplace(chrom, (uint32_t)chrom_index.size());
        if (result.second) {
            chrom_is_gref_fragment.push_back(GrefCover::is_gref_name(chrom));
        }
        return result.first->second;
    };
    // One entry per snarl name.  A snarl ID is not unique -- a cyclic reference path that
    // traverses the same snarl twice emits two records with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa) -- but both occurrences are traversals of one
    // path, so they are on the same contig and it does not matter which one wins here.
    unordered_map<string, uint32_t> chrom_of_name;
    // The site name and reference stretch of each record, decompressed once, for the indexes and
    // the second pass.
    vector<vector<pair<string, Region>>> placed = read_record_regions(output_variants);
    // Fill the indexes in record order.
    for (const vector<pair<string, Region>>& buffer : placed) {
        for (const auto& [site_name, region] : buffer) {
            chrom_of_name.emplace(site_name, intern_chrom(region.seq));
        }
    }

    // Index the sites by the names of their records. Only names in the VCF are looked up, so only
    // they are indexed.
    unordered_map<string, net_handle_t> name_to_snarl =
        index_sites_by_record_name(sites, chrom_of_name, record_id_of_site);

    // pass 2) identify top-level snarls (those with no ancestors in VCF)
    // and store reference info only for them
    // Keyed by snarl name, then by (chrom, pos), because a snarl ID can carry more than one
    // record: a cyclic reference emits two, both with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa, which gives two <5<1 records at POS 20 and 44).
    // A plain name -> interval map was last-write-wins, so both records were handed the
    // surviving one's interval and the record at POS 20 reported RS=44.  The inner map is
    // ordered so that picking begin() is deterministic regardless of thread scheduling.
    //
    // The inner map gives the end of each record's reference stretch by its contig and start.
    unordered_map<string, map<pair<string, int64_t>, int64_t>> top_level_ref_info;

    // Helper to check if a snarl is top-level (no ancestors in VCF)
    auto is_top_level = [&](const string& name) {
        auto it = name_to_snarl.find(name);
        if (it == name_to_snarl.end()) return true; // not found, treat as top-level
        for (net_handle_t snarl = parent_site(sites, it->second); !sites.is_root(snarl);
             snarl = parent_site(sites, snarl)) {
            string cur_name = name_of(snarl);
            string flipped_name = flipped_name_of(snarl);
            if (chrom_of_name.contains(cur_name) || chrom_of_name.contains(flipped_name)) {
                return false; // has ancestor in VCF
            }
        }
        return true; // no ancestors in VCF
    };

    // Find the records whose sites are top level, in parallel. Then store their reference
    // intervals in record order.
    // Whether the site of each record is top level, by buffer.
    vector<vector<bool>> top_level(placed.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t b = 0; b < placed.size(); ++b) {
        for (const auto& [site_name, region] : placed[b]) {
            top_level[b].push_back(is_top_level(site_name));
        }
    }
    for (size_t b = 0; b < placed.size(); ++b) {
        for (size_t i = 0; i < placed[b].size(); ++i) {
            if (top_level[b][i]) {
                const auto& [site_name, region] = placed[b][i];
                top_level_ref_info[site_name][make_pair(region.seq, region.start)] = region.end;
            }
        }
    }
    vector<vector<pair<string, Region>>>().swap(placed);

    // determine the tags from the index
    //
    // There are exactly two ways a snarl can nest inside its parent's record, and they need
    // to be told apart.  A site inside a *deletion* is covered by its parent contig's own
    // reference allele, so it has coordinates on that contig and its record's CHROM is the
    // same.  A site inside an *insertion* has no path of the parent's contig through it at
    // all, so it is only callable once some other reference (a gref fragment) covers the
    // inserted allele -- and its record's CHROM is therefore different.  So:
    //
    //   contig_level    ancestors whose record is on this record's own CHROM, i.e. how deep
    //                   the site is in its own coordinate system
    //   contig_hops     steps in the chain where CHROM changed, i.e. how many insertions deep
    //                   the site is
    //
    // Returns: (contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
    //           suppressed_name)
    // ref_chrom_name is the topmost ancestor in the VCF that sits on a reference contig rather
    // than a gref one, and suppressed_name the topmost ancestor that was dropped for having no
    // variant.  Both feed the RC/RS/RD choice below; neither affects LV/CH/PS.
    function<tuple<size_t, size_t, string, string, string, string>(const string&, const string&)> get_nesting_tags =
        [&](const string& name, const string& my_chrom) {
        string parent_name;
        string ref_chrom_name;
        string suppressed_name;
        string top_level_name = name;  // default to self (for the top-level case)
        size_t contig_level = 0;
        size_t contig_hops = 0;
        // Our own contig, and the contig of the previously visited link in the chain.
        // chrom_index is complete after pass 1, so this lookup always hits.
        uint32_t my_chrom_id = chrom_index.at(my_chrom);
        uint32_t prev_chrom_id = my_chrom_id;
        // walk up the snarl tree
        for (net_handle_t snarl = parent_site(sites, name_to_snarl.at(name));
             !sites.is_root(snarl); snarl = parent_site(sites, snarl)) {
            string cur_name = name_of(snarl);

            // Since it is possible that the snarl is actually flipped in the vcf, check for the
            // flipped version too
            string flipped_name = flipped_name_of(snarl);
            const string* hit = nullptr;
            if (chrom_of_name.contains(cur_name)) {
                // only count snarls that are in the vcf
                hit = &cur_name;
            } else if (chrom_of_name.contains(flipped_name)) {
                // snarl is in vcf under flipped orientation
                hit = &flipped_name;
            }
            if (hit == nullptr) {
                // Not in the VCF.  If it was dropped for having no variant we still know where it
                // sits, and it may be the only ancestor that can give a reference position.
                auto sup_it = suppressed_ref.find(cur_name);
                if (sup_it == suppressed_ref.end()) {
                    sup_it = suppressed_ref.find(flipped_name);
                }
                if (sup_it != suppressed_ref.end()) {
                    suppressed_name = sup_it->first;
                }
                continue;
            }

            auto chrom_it = chrom_of_name.find(*hit);
            uint32_t anc_chrom_id = chrom_it == chrom_of_name.end() ? my_chrom_id
                                                                    : chrom_it->second;
            if (anc_chrom_id == my_chrom_id) {
                ++contig_level;
            }
            if (anc_chrom_id != prev_chrom_id) {
                ++contig_hops;
            }
            prev_chrom_id = anc_chrom_id;

            if (parent_name.empty()) {
                // remember the first parent
                parent_name = *hit;
            }
            // keep updating top_level to find the topmost ancestor in VCF
            top_level_name = *hit;
            // ...and, separately, the topmost one actually on a reference contig.  An ancestor on
            // another gref contig can name a position, but not one a reader can look up in the
            // reference, so it is the weaker answer of the two.
            if (!chrom_is_gref_fragment[anc_chrom_id]) {
                ref_chrom_name = *hit;
            }
        }
        return make_tuple(contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                          suppressed_name);
    };

    // pass 3) add the LV, PS, RC, RS, RD tags
#pragma omp parallel for
    for (uint64_t i = 0; i < output_variants.size(); ++i) {
        auto& thread_buf = output_variants[i];
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            // Look the record up by the name of its site, so that a block record gets the tags of
            // its site.
            const string name = block_site_name(toks[2]);

            auto [contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                  suppressed_name] = get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig, so that a gRef fragment's
            // records start at level 0 on their own contig.
            //
            // CH counts the ancestors with a record, which is 0 on a gRef fragment whose enclosing
            // site has no line. So we use the gRef level of the contig as a floor for CH. Then
            // CH >= 1 does not imply a parent record, or PS.
            size_t gref_level = 0;
            {
                auto it = gref_levels.find(toks[0]);
                if (it != gref_levels.end() && it->second > 0) {
                    gref_level = (size_t)it->second;
                }
            }
            string nesting_tags = ";LV=" + std::to_string(contig_level);
            nesting_tags += ";CH=" + std::to_string(max(contig_hops, gref_level));
            if (!parent_name.empty()) {
                // Not "if (lv != 0)": those were equivalent only while LV was the absolute
                // count.  A record can now legitimately be at LV=0 and still have a parent on
                // another contig, and it must keep PS -- vcfbub's rescue of the children of
                // popped bubbles is keyed on it.
                nesting_tags += ";PS=" + parent_name;
            }

            // Add RC, RS, RD tags: where to look this record up in the reference.
            //
            // Prefer, in order, the topmost ancestor in the VCF that is on a reference contig;
            // then the topmost ancestor in the VCF at all; then the topmost ancestor that was
            // dropped for having no variant.  The last is what rescues a gref fragment whose
            // parent snarl only the reference and its own gref copy span: the site is real and
            // has a reference interval, it just had nothing to report.
            //
            // If none of those exist there is no reference position to give, and the tags are
            // left off.  They used to fall back to this record's own contig and position, which
            // is not a reference coordinate at all -- on a gref contig it is a self-reference
            // that a reader cannot tell apart from the genuine case.
            const string* ref_source = nullptr;
            if (!ref_chrom_name.empty()) {
                ref_source = &ref_chrom_name;
            } else if (top_level_name != name) {
                ref_source = &top_level_name;
            }
            bool have_ref = true;
            // The stretch of the reference to report, 1-based and end-inclusive.
            Region top_ref;
            if (ref_source == nullptr) {
                if (!GrefCover::is_gref_name(toks[0])) {
                    // Not on a gref fragment, so our own interval is already a position a reader
                    // can look up, and it is the narrower answer of the two.  Keep it rather than
                    // reach for an enclosing site that produced no record: doing that would
                    // repoint every such record on a reference contig at a site LV and CH say it
                    // has no ancestor in.  The records that need the reach are the fragments
                    // below, which have no usable coordinate of their own.
                    int64_t position = stoll(toks[1]);
                    top_ref = {toks[0], position, position + (int64_t)toks[3].length() - 1};
                } else {
                    auto sup_it = suppressed_name.empty() ? suppressed_ref.end()
                                                          : suppressed_ref.find(suppressed_name);
                    if (sup_it != suppressed_ref.end()) {
                        top_ref = sup_it->second;
                    } else {
                        have_ref = false;
                    }
                }
            } else {
                const auto& candidates = top_level_ref_info.at(*ref_source);
                // If the ancestor produced several records, prefer one on our own contig;
                // failing that take the smallest (chrom, pos).  Which one is "right" is
                // genuinely ambiguous, so pick deterministically rather than by chance.
                auto chosen = candidates.begin();
                for (auto it = candidates.begin(); it != candidates.end(); ++it) {
                    if (it->first.first == toks[0]) {
                        chosen = it;
                        break;
                    }
                }
                top_ref = {chosen->first.first, chosen->first.second, chosen->second};
            }
            if (have_ref) {
                nesting_tags += ";RC=" + top_ref.seq;
                nesting_tags += ";RS=" + std::to_string(top_ref.start);
                nesting_tags += ";RD=" + std::to_string(top_ref.end + 1);
            }

            // rewrite the output string using the updated info toks
            output_variant_string.clear();
            for (size_t i = 0; i < toks.size(); ++i) {
                output_variant_string += toks[i];
                if (i == 7) {
                    output_variant_string += nesting_tags;
                }
                if (i != toks.size() - 1) {
                    output_variant_string += "\t";
                }
            }
            output_variant_record.second.clear();
            ret = zstdutil::CompressString(output_variant_string, output_variant_record.second);
            assert(ret == 0);
        }
    }
}
}

