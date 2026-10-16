#include "vcf_output_caller.hpp"
#include "algorithms/expand_context.hpp"
#include "annotation.hpp"
#include "gref.hpp"
#include "traversal_clusters.hpp"

//#define debug

namespace vg {
    
VCFOutputCaller::VCFOutputCaller(const string& sample_name) : sample_name(sample_name), translation(nullptr), include_nested(false)
{
    output_variants.resize(get_thread_count());
    suppressed_ref_info.resize(get_thread_count());
}

VCFOutputCaller::~VCFOutputCaller() {
}

string VCFOutputCaller::vcf_header(const PathHandleGraph& graph, const vector<string>& contigs,
                                   const vector<size_t>& contig_length_overrides) const {
    stringstream ss;
    ss << "##fileformat=VCFv4.2" << endl;    
    for (int i = 0; i < contigs.size(); ++i) {
        const string& contig = contigs[i];
        size_t length;
        if (i < contig_length_overrides.size()) {
            // length override provided
            length = contig_length_overrides[i];
        } else {
            length = 0;
            for (handle_t handle : graph.scan_path(graph.get_path_handle(contig))) {
                length += graph.get_length(handle);
            }
        }
        ss << "##contig=<ID=" << contig << ",length=" << length << ">" << endl;
    }
    if (include_nested) {
        ss << nesting_info_headers();
    }
    ss << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele Traversal as path in graph\">" << endl;
    if (allele_merge_threshold < 1.0) {
        ss << "##INFO=<ID=MAT,Number=.,Type=String,Description=\"Merged Allele Traversal: "
           << "ALT alleles merged after genotyping by -L/--cluster, as OLD>NEW:SIMILARITY using "
           << "pre-merge allele numbers. In a nested run this record gives the collapsed view of "
           << "the site and its child records the precise one, so they disagree by design. "
           << "pre-merge allele numbers. AD and GL are folded onto the surviving allele and MAD is "
           << "recomputed; DP, QUAL, GQ, GP and FILTER are as computed over the pre-merge allele set.\">"
           << endl;
    }
    return ss.str();
}

bool VCFOutputCaller::add_variant(vcflib::Variant& var) const {
    var.setVariantCallFile(output_vcf);
    stringstream ss;
    ss << var;
    string dest;
    if (ss.str().length() > VCFOutputCaller::max_vcf_line_length) {
        return false;
    }         
    int ret = zstdutil::CompressString(ss.str(), dest);
    assert(ret == 0);
    // the Variant object is too big to keep in memory when there are many genotypes, so we
    // store it in a zstd-compressed string
    output_variants[omp_get_thread_num()].push_back(make_pair(make_pair(var.sequenceName, var.position), dest));
    return true;
}

void VCFOutputCaller::write_variants(ostream& out_stream, const SnarlManager* snarl_manager) {
    assert(include_nested == false || snarl_manager != nullptr);
    if (include_nested) {
        update_nesting_info_tags(snarl_manager);
    }
    vector<pair<pair<string, size_t>, string>> all_variants;
    // Reserve once: doing it inside the loop below reallocates per thread buffer.
    size_t total_variants = 0;
    for (const auto& buf : output_variants) {
        total_variants += buf.size();
    }
    all_variants.reserve(total_variants);
    // `buf` must not be const: std::move() over const_iterators silently degrades to a copy,
    // which duplicated every compressed record at the one point where the whole VCF is in
    // memory at once.  Free the buffer as we go for the same reason.
    //
    // This makes write_variants() single-use, which it already effectively was -- a real move
    // leaves the buffers empty either way.  All three callers (deconstructor.cpp,
    // call_main.cpp, mcmc_main.cpp) call it exactly once.
    for (auto& buf : output_variants) {
        std::move(buf.begin(), buf.end(), std::back_inserter(all_variants));
        buf.clear();
        buf.shrink_to_fit();
    }
    std::sort(all_variants.begin(), all_variants.end(), [](const pair<pair<string, size_t>, string>& v1,
                                                           const pair<pair<string, size_t>, string>& v2) {
            return v1.first.first < v2.first.first || (v1.first.first == v2.first.first && v1.first.second < v2.first.second);
        });
    for (const auto& v : all_variants) {
        string dest;
        int ret = zstdutil::DecompressString(v.second, dest);
        assert(ret == 0);
        out_stream << dest << endl;
    }
}

static int countAlts(vcflib::Variant& var, int alleleIndex) {
    int alts = 0;
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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
    for (map<string, map<string, vector<string> > >::iterator s = var.samples.begin(); s != var.samples.end(); ++s) {
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

void VCFOutputCaller::vcf_fixup(vcflib::Variant& var) const {
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

void VCFOutputCaller::set_translation(const unordered_map<nid_t, pair<string, size_t>>* translation) {
    this->translation = translation;
}

void VCFOutputCaller::set_nested(bool nested) {
    include_nested = nested;
}

void VCFOutputCaller::set_allele_merge(double threshold, int64_t min_len) {
    allele_merge_threshold = threshold;
    allele_merge_min_len = min_len;
}

bool VCFOutputCaller::snarl_traversal_to_handles(const HandleGraph& graph, const SnarlTraversal& trav,
                                                 Traversal& out_trav) {
    // cluster_traversals asserts size() >= 2, and a Visit carrying a child Snarl has no single
    // handle.  Both are real inputs here (the "*" placeholder, and NestedFlowCaller traversals via
    // SnarlGraph::embed_snarl), so refuse rather than fabricate something.
    if (trav.visit_size() < 2) {
        return false;
    }
    out_trav.clear();
    out_trav.reserve(trav.visit_size());
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() <= 0) {
            return false;
        }
        out_trav.push_back(graph.get_handle(visit.node_id(), visit.backward()));
    }
    return true;
}

namespace {
/// Parse a VCF FORMAT value without throwing and without exiting.  vg::parse<double> exits on
/// failure and the 2-argument vg::parse can throw; merge_similar_alleles runs inside an OpenMP
/// region, where an escaping exception is std::terminate rather than something a caller can handle,
/// and exit() would abandon whatever the other threads had already buffered.  Missing values
/// ("." and "") are ordinary input here.
bool parse_vcf_double(const string& field, double& value) {
    try {
        size_t after;
        value = std::stod(field, &after);
        return after == field.size();
    } catch (const std::exception&) {
        return false;
    }
}
}

int64_t VCFOutputCaller::allele_core_length(const vector<string>& alleles) {
    vector<const string*> seqs;
    for (const string& a : alleles) {
        if (a != "*") {
            seqs.push_back(&a);
        }
    }
    if (seqs.empty()) {
        return 0;
    }
    size_t min_len = seqs[0]->length();
    size_t max_len = 0;
    for (const string* s : seqs) {
        min_len = std::min(min_len, s->length());
        max_len = std::max(max_len, s->length());
    }
    // The prefix and the suffix may not overlap, exactly as in flatten_common_allele_ends: a shared
    // region can only be counted once, or {"AAAA","AAAAA"} would come out at -3 instead of 1.
    // Case-insensitive to match flatten's own toupper and deconstruct's toUppercase.
    auto shared = [&](size_t skip, bool from_back) {
        auto at = [&](const string* s, size_t i) {
            return std::toupper((*s)[from_back ? s->length() - 1 - i : i]);
        };
        size_t n = 0;
        while (skip + n < min_len) {
            int ch = at(seqs[0], n);
            bool match = true;
            for (size_t j = 1; j < seqs.size() && match; ++j) {
                match = at(seqs[j], n) == ch;
            }
            if (!match) {
                break;
            }
            ++n;
        }
        return n;
    };
    size_t prefix = shared(0, false);
    size_t suffix = shared(prefix, true);
    // non-negative structurally, not by clamping: the loop caps give prefix + suffix <= min_len <= max_len
    return (int64_t)(max_len - prefix - suffix);
}

bool VCFOutputCaller::merge_similar_alleles(const PathPositionHandleGraph& graph,
                                            const vector<SnarlTraversal>& site_traversals,
                                            vector<int>& site_genotype,
                                            const string& sample_name,
                                            vcflib::Variant& out_variant) const {
    if (!(allele_merge_threshold < 1.0)) {
        return false;
    }
    // we only collapse a genotype that actually calls two distinct ALTs.  This also keeps the -a
    // padding block above (which adds uncalled alleles) out of scope: those alleles are advertised,
    // not called, and rewriting them without a genotype change would be a silent surprise.
    set<int> called_alts;
    for (int g : site_genotype) {
        if (g > 0) {
            called_alts.insert(g);
        }
    }
    if (called_alts.size() < 2) {
        return false;
    }
    // Per-site gate, decided over the alleles this record actually emits.  NOT over the traversal
    // finder's candidate list: that is up to max_yens_traversals (50) speculative paths, most of
    // which never become an allele, so gating on them lets an invisible branch with no reads and no
    // AT entry decide whether merging happens.  deconstruct's equivalent gate is decided over the
    // set that becomes ITS alleles -- the reference plus everything get_traversal_order kept.  Both
    // tools gate on what they emit, but those sets differ (we see only the called genotype's
    // alleles, deconstruct sees every haplotype), so the two can disagree at a site whose uncalled
    // haplotypes are much larger than its called ones.
    // The quantity is CORE LENGTH (see allele_core_length): the longest allele once the prefix and
    // suffix shared by every allele are stripped.  Raw string length would answer differently from
    // vg deconstruct on the same variant, because this record has been flattened down to an anchor
    // base and deconstruct's has not.
    if (allele_merge_min_len > 0 &&
        allele_core_length(out_variant.alleles) < allele_merge_min_len) {
        return false;
    }

    // ALT-vs-ALT only.  Absorbing an ALT into allele 0 would empty out_variant.alt and the record
    // would then be dropped entirely by the caller, turning a het call into no call at all.
    // (vg deconstruct does fold near-reference alleles into the reference cluster and drop the
    //  record; that is deliberate there and deliberately not copied here.)
    vector<Traversal> alt_travs;
    vector<int> alt_to_allele;
    for (size_t i = 1; i < site_traversals.size(); ++i) {
        if (!called_alts.count((int)i)) {
            continue;
        }
        Traversal trav;
        if (!snarl_traversal_to_handles(graph, site_traversals[i], trav)) {
            // star placeholder or a child-snarl visit: leave this allele alone
            continue;
        }
        alt_travs.push_back(std::move(trav));
        alt_to_allele.push_back((int)i);
    }
    if (alt_travs.size() < 2) {
        return false;
    }

    // same clustering call deconstruct makes, so the metric, the endpoint pruning and the
    // >= comparison are inherited rather than reimplemented.
    //
    // Cluster in descending allele-depth order, so each cluster's head -- the allele that survives,
    // and the one MAT's similarity is measured against -- is its best-supported member.  Identity order
    // would instead inherit the traversal finder's ranking, which FlowCaller::call_snarl_internal (and NestedFlowCaller's copy of it) switches to
    // length-weighted average flow once a snarl's interior passes the average-support threshold.
    // That ranking can put a short, lightly-supported allele ahead of a long, heavily-supported one,
    // and merging into it emits the minority sequence as a homozygous call carrying the pooled depth.
    vector<int> order(alt_travs.size());
    std::iota(order.begin(), order.end(), 0);
    {
        auto& sample_fields = out_variant.samples[sample_name];
        auto ad_it = sample_fields.find("AD");
        if (ad_it != sample_fields.end() && ad_it->second.size() == out_variant.alleles.size()) {
            vector<double> ad(alt_travs.size(), 0);
            bool usable = true;
            for (size_t k = 0; k < alt_travs.size() && usable; ++k) {
                usable = parse_vcf_double(ad_it->second.at(alt_to_allele[k]), ad[k]);
            }
            if (usable) {
                // stable, so equal depths keep the finder's own ranking
                std::stable_sort(order.begin(), order.end(),
                                 [&](int a, int b) { return ad[a] > ad[b]; });
            }
        }
    }
    // The VCF reference allele sets the scale a pure deletion is measured against.  It is
    // site_traversals[0] and is deliberately not clustered (the loop above starts at 1).
    // The nullptr fallback is currently unreachable: the only producer of a Visit without a node is
    // NestedFlowCaller, and --bottom-up is rejected with -L.  It is kept because the consequence of
    // being wrong about that is a crash, and because falling back to pairwise scoring can only
    // merge less, never more.
    Traversal ref_trav;
    const Traversal* site_ref_trav = nullptr;
    if (!site_traversals.empty() && snarl_traversal_to_handles(graph, site_traversals[0], ref_trav)) {
        site_ref_trav = &ref_trav;
    }
    vector<pair<double, int64_t>> cluster_info;
    vector<int> unused_child_snarl_mapping;
    vector<vector<int>> clusters = cluster_traversals(&graph, alt_travs, order,
                                                      vector<pair<handle_t, handle_t>>(),
                                                      allele_merge_threshold,
                                                      cluster_info, unused_child_snarl_mapping,
                                                      site_ref_trav);

    // merge_to[a] == a for a surviving allele, else the allele it collapses into
    vector<int> merge_to(out_variant.alleles.size());
    std::iota(merge_to.begin(), merge_to.end(), 0);
    vector<string> mat_entries;
    bool merged_any = false;
    for (const vector<int>& cluster : clusters) {
        if (cluster.size() < 2) {
            continue;
        }
        int survivor = alt_to_allele[cluster.front()];
        for (size_t j = 1; j < cluster.size(); ++j) {
            int absorbed = alt_to_allele[cluster[j]];
            merge_to[absorbed] = survivor;
            merged_any = true;
            stringstream ss;
            ss.precision(3);
            ss << absorbed << ">" << survivor << ":" << cluster_info[cluster[j]].first;
            mat_entries.push_back(ss.str());
        }
    }
    if (!merged_any) {
        return false;
    }

    // dense renumbering of the survivors, preserving order
    vector<int> new_index(merge_to.size(), -1);
    int next = 0;
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_index[a] = next++;
        }
    }
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] != (int)a) {
            new_index[a] = new_index[merge_to[a]];
        }
    }
    int n_new = next;

    // alleles / alt
    vector<string> new_alleles(n_new);
    for (size_t a = 0; a < merge_to.size(); ++a) {
        if (merge_to[a] == (int)a) {
            new_alleles[new_index[a]] = out_variant.alleles[a];
        }
    }
    out_variant.alleles = new_alleles;
    out_variant.alt.assign(new_alleles.begin() + 1, new_alleles.end());

    // AT is Number=R, so it is indexed by allele just like alleles
    auto at_it = out_variant.info.find("AT");
    if (at_it != out_variant.info.end() && at_it->second.size() == merge_to.size()) {
        vector<string> new_at(n_new);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            if (merge_to[a] == (int)a) {
                new_at[new_index[a]] = at_it->second[a];
            }
        }
        at_it->second = new_at;
    }

    auto& sample = out_variant.samples[sample_name];

    // AD is a per-allele count, so the absorbed allele's reads move onto the survivor.  sum(AD) is
    // therefore unchanged by the merge; it can slightly over-count when the merged alleles share
    // interior nodes, whose depth was proportionally split between them.
    auto ad_it = sample.find("AD");
    if (ad_it != sample.end() && ad_it->second.size() == merge_to.size()) {
        vector<double> summed(n_new, 0);
        for (size_t a = 0; a < merge_to.size(); ++a) {
            double v = 0;
            // treat an unparseable entry as 0 rather than bailing: the merge is already committed,
            // and dropping AD would leave the record with no per-allele depth at all
            parse_vcf_double(ad_it->second[a], v);
            summed[new_index[a]] += v;
        }
        vector<string> new_ad(n_new);
        for (int a = 0; a < n_new; ++a) {
            new_ad[a] = std::to_string((int64_t)std::llround(summed[a]));
        }
        ad_it->second = new_ad;
        // MAD is the min allele depth over the called alleles; recompute so it agrees with the AD
        // and GT printed beside it.  FILTER is left as computed pre-merge, and since the new MAD is
        // >= the old one, a depth filter can only over-filter, never under-filter.
        auto mad_it = sample.find("MAD");
        if (mad_it != sample.end() && mad_it->second.size() == 1) {
            double min_ad = -1;
            for (int g : site_genotype) {
                if (g >= 0 && g < (int)merge_to.size()) {
                    double v = summed[new_index[g]];
                    if (min_ad < 0 || v < min_ad) {
                        min_ad = v;
                    }
                }
            }
            if (min_ad >= 0) {
                mad_it->second[0] = std::to_string((int64_t)std::llround(min_ad));
            }
        }
    }

    // GL is Number=G.  vg emits it i-major -- "for i; for j = i..n" -- which is not the VCF spec's
    // ordering for 3+ alleles, but the fold has to match what is actually written.  Take the max
    // over the old genotype classes mapping onto each new one: that is the max-marginal, i.e. the
    // merged allele scores as whichever of its members fit best.
    auto gl_it = sample.find("GL");
    if (gl_it != sample.end()) {
        size_t n_old = merge_to.size();
        auto gl_index = [](size_t i, size_t j, size_t n) { return i * n - (i * (i - 1)) / 2 + (j - i); };
        // The diploid layout is the only one that can occur: merging needs at least two distinct
        // called ALTs, site_genotype has one entry per ploidy, and PoissonSupportSnarlCaller --
        // whose update_vcf_info is the only writer of GL -- asserts in genotype that ploidy is
        // 1 or 2.  So n_old is always 3 here.
        assert(gl_it->second.size() == n_old * (n_old + 1) / 2);
        bool gl_usable = true;
        vector<double> folded((size_t)n_new * (n_new + 1) / 2,
                              -std::numeric_limits<double>::infinity());
        for (size_t i = 0; i < n_old && gl_usable; ++i) {
            for (size_t j = i; j < n_old && gl_usable; ++j) {
                double v = 0;
                gl_usable = parse_vcf_double(gl_it->second[gl_index(i, j, n_old)], v);
                if (!gl_usable) {
                    break;
                }
                size_t ni = new_index[i], nj = new_index[j];
                if (ni > nj) {
                    std::swap(ni, nj);
                }
                double& slot = folded[gl_index(ni, nj, (size_t)n_new)];
                slot = std::max(slot, v);
            }
        }
        if (gl_usable) {
            vector<string> new_gl(folded.size());
            for (size_t i = 0; i < folded.size(); ++i) {
                new_gl[i] = std::to_string(folded[i]);
            }
            gl_it->second = new_gl;
        } else {
            // A value we cannot parse.  Leaving GL alone would emit a Number=G field whose length
            // disagrees with the new allele count, so drop it rather than lie.
            sample.erase(gl_it);
            auto& fmt = out_variant.format;
            fmt.erase(std::remove(fmt.begin(), fmt.end(), string("GL")), fmt.end());
        }
    }
    // GQ and GP are deliberately untouched: they come from the caller's own CallInfo, computed over
    // its candidate set rather than from the emitted GL, so recomputing them here would silently
    // swap one statistic for another.

    // GT, from the renumbered genotype
    for (int& g : site_genotype) {
        if (g >= 0 && g < (int)merge_to.size()) {
            g = new_index[g];
        }
    }
    stringstream vcf_gt;
    for (size_t i = 0; i < site_genotype.size(); ++i) {
        if (site_genotype[i] == MISSING_ALLELE_MARKER) {
            vcf_gt << ".";
        } else {
            vcf_gt << site_genotype[i];
        }
        if (i != site_genotype.size() - 1) {
            vcf_gt << "/";
        }
    }
    sample["GT"] = {vcf_gt.str()};

    // record what happened: without this a merged 1/1 is indistinguishable from a real hom-alt
    out_variant.info["MAT"] = mat_entries;

    out_variant.updateAlleleIndexes();
    return true;
}

unordered_set<string> VCFOutputCaller::get_output_contigs() const {
    unordered_set<string> contigs;
    // The sort key is (sequenceName, position) (see add_variant), so the contig is right
    // there and nothing has to be decompressed.
    for (const auto& thread_buf : output_variants) {
        for (const auto& output_variant_record : thread_buf) {
            contigs.insert(output_variant_record.first.first);
        }
    }
    return contigs;
}

string VCFOutputCaller::prune_header_contigs(const string& header,
                                             const unordered_set<string>& keep) const {
    static const string contig_prefix = "##contig=<ID=";
    stringstream pruned;
    vector<string> lines = split_delims(header, "\n");
    for (const string& line : lines) {
        if (line.compare(0, contig_prefix.size(), contig_prefix) == 0) {
            // Parse the ID back out the same way it was written, rather than scanning for a
            // delimiter: contig names are path names and nothing stops one containing ',' or
            // '>'.  Both producers emit exactly ##contig=<ID=NAME,length=N> -- vcf_header()
            // above and Deconstructor::add_contigs_to_vcf_header().
            static const string contig_suffix = ",length=";
            size_t id_start = contig_prefix.size();
            size_t id_end = line.rfind(contig_suffix);
            if (id_end == string::npos || id_end < id_start) {
                // not a shape we wrote; leave it alone rather than guess
                pruned << line << "\n";
                continue;
            }
            string id = line.substr(id_start, id_end - id_start);
            if (!keep.count(id)) {
                continue;
            }
        }
        pruned << line << "\n";
    }
    string result = pruned.str();
    if (!header.empty() && header.back() != '\n' && !result.empty()) {
        // input had no trailing newline, so don't invent one
        result.pop_back();
    }
    return result;
}

void VCFOutputCaller::add_allele_path_to_info(const HandleGraph* graph, vcflib::Variant& v, int allele, const Traversal& trav,
                                              bool reversed, bool one_based) const {
    SnarlTraversal proto_trav;
    for (const handle_t& handle : trav) {
        Visit* visit = proto_trav.add_visit();
        visit->set_node_id(graph->get_id(handle));
        visit->set_backward(graph->get_is_reverse(handle));
    }
    this->add_allele_path_to_info(v, allele, proto_trav, reversed, one_based);
}

void VCFOutputCaller::add_allele_path_to_info(vcflib::Variant& v, int allele, const SnarlTraversal& trav,
                                              bool reversed, bool one_based) const {
    auto& trav_info = v.info["AT"];
    assert(allele < trav_info.size());

    vector<int> nodes;
    nodes.reserve(trav.visit_size());
    const Visit* prev_visit = nullptr;
    unordered_map<nid_t, pair<string, size_t>>::const_iterator prev_trans;
    
    for (size_t i = 0; i < trav.visit_size(); ++i) {
        size_t j = !reversed ? i : trav.visit_size() - 1 - i;
        const Visit& visit = trav.visit(j);
        nid_t node_id = visit.node_id();
        string node_name = std::to_string(node_id);
        bool skip = false;
        // todo: check one_based? (we kind of ignore that when writing the snarl name, so maybe not pertienent)
        if (translation) {
            auto i = translation->find(node_id);
            if (i == translation->end()) {
                throw runtime_error("Error [vg deconstruct]: Unable to find node " + node_name + " in translation file");
            }
            if (prev_visit) {
                nid_t prev_node_id = prev_visit->node_id();
                if (prev_trans->second.first == i->second.first && node_id != prev_node_id) {
                    // here is a case where we have two consecutive nodes that map back to
                    // the same source node.
                    // todo: check if translation node properly covered
                    skip = true;
                }
            }
            node_name = i->second.first;
            prev_trans = i;
        }

        if (!skip) {
            bool vrev = visit.backward() != reversed;
            trav_info[allele] += (vrev ? "<" : ">");
            trav_info[allele] += node_name;
        }
        prev_visit = &visit;
    }
    if (trav_info[allele].empty()) {
        // note: * alleles get empty traversals
        trav_info[allele] = ".";
    }
}

string VCFOutputCaller::trav_string(const HandleGraph& graph, const SnarlTraversal& trav) const {
    string seq;
    for (int i = 0; i < trav.visit_size(); ++i) {
        const Visit& visit = trav.visit(i);
        if (visit.node_id() > 0) {
            seq += graph.get_sequence(graph.get_handle(visit.node_id(), visit.backward()));
        } else {
            seq += print_snarl(visit.snarl(), true);
        }
    }
    return seq;    
}

bool VCFOutputCaller::emit_variant(const PathPositionHandleGraph& graph, SnarlCaller& snarl_caller,
                                   const Snarl& snarl, const vector<SnarlTraversal>& called_traversals,
                                   const vector<int>& genotype, int ref_trav_idx, const unique_ptr<SnarlCaller::CallInfo>& call_info,
                                   const string& ref_path_name, int ref_offset, bool genotype_snarls, int ploidy,
                                   function<string(const vector<SnarlTraversal>&, const vector<int>&, int, int, int)> trav_to_string) {
    
#ifdef debug
    cerr << "emitting variant for " << pb2json(snarl) << endl;
    for (int i = 0; i < called_traversals.size(); ++i) {
        if (i == ref_trav_idx) {
            cerr << "*";
        }
        cerr << "ct[" << i << "]=" << pb2json(called_traversals[i]) << endl;
    }
    for (int i = 0; i < genotype.size(); ++i) {
        cerr << "gt[" << i << "]=" << genotype[i] << endl;
    }
#endif

    if (trav_to_string == nullptr) {
        trav_to_string = [&](const vector<SnarlTraversal>& travs, const vector<int>& travs_genotype, int trav_allele, int genotype_allele, int ref_trav_idx) {
            return trav_string(graph, travs[trav_allele]);    
        };
    }

    vcflib::Variant out_variant;

    vector<SnarlTraversal> site_traversals = {called_traversals[ref_trav_idx]};
    vector<int> site_genotype;
    auto ref_gt_it = std::find(genotype.begin(), genotype.end(), ref_trav_idx);
    out_variant.ref = trav_to_string(called_traversals, genotype, ref_trav_idx,
                                     ref_gt_it != genotype.end() ? ref_gt_it - genotype.begin() : 0,
                                     ref_trav_idx);
    
    // deduplicate alleles and compute the site traversals and genotype
    map<string, int> allele_to_gt;
    allele_to_gt[out_variant.ref] = 0;
    int star_allele_idx = -1;  // index for star allele in allele_to_gt, if needed
    for (int i = 0; i < genotype.size(); ++i) {
        if (genotype[i] == STAR_ALLELE_MARKER) {
            // Star allele: haplotype spans this site but has no defined traversal here
            if (star_allele_idx < 0) {
                // Add star allele to allele list
                star_allele_idx = allele_to_gt.size();
                allele_to_gt["*"] = star_allele_idx;
                // Add empty traversal as placeholder (won't be used for AT info)
                site_traversals.push_back(SnarlTraversal());
            }
            site_genotype.push_back(star_allele_idx);
        } else if (genotype[i] == MISSING_ALLELE_MARKER) {
            // Missing allele: parent doesn't traverse this child, output as '.' in VCF
            site_genotype.push_back(MISSING_ALLELE_MARKER);
        } else if (genotype[i] == ref_trav_idx) {
            site_genotype.push_back(0);
        } else {
            string allele_string = trav_to_string(called_traversals, genotype, genotype[i], i, ref_trav_idx);
            if (allele_to_gt.count(allele_string)) {
                site_genotype.push_back(allele_to_gt[allele_string]);
            } else {
                site_traversals.push_back(called_traversals[genotype[i]]);
                site_genotype.push_back(allele_to_gt.size());
                allele_to_gt[allele_string] = site_genotype.back();
            }
        }
    }

    // add on fixed number of uncalled traversals if we're making a ref-call
    // with genotype_snarls set to true
    if (genotype_snarls && site_traversals.size() <= 1) {
        // note: we're adding all the strings here and sorting to make this deterministic
        // at the cost of speed
        map<string, const SnarlTraversal*> allele_map;
        for (int i = 0; i < called_traversals.size(); ++i) {
            // todo: verify index below.  it's for uncalled traversals so not important tho
            string allele_string = trav_to_string(called_traversals, genotype, i, max(0, (int)genotype.size() - 1), ref_trav_idx);
            if (!allele_map.count(allele_string)) {
                allele_map[allele_string] = &called_traversals[i];
            }
        }
        // pick out the first "max_uncalled_alleles" traversals to add
        int i = 0;
        for (auto ai = allele_map.begin(); i < max_uncalled_alleles && ai != allele_map.end(); ++i, ++ai) {
            if (!allele_to_gt.count(ai->first)) {
                allele_to_gt[ai->first] = allele_to_gt.size();
                site_traversals.push_back(*ai->second);
            }
        }
    }

    out_variant.alt.resize(allele_to_gt.size() - 1);
    out_variant.alleles.resize(allele_to_gt.size());
    
    // init the traversal info
    out_variant.info["AT"].resize(allele_to_gt.size());

    for (auto& allele_gt : allele_to_gt) {
#ifdef debug
        cerr << "allele " << allele_gt.first << " -> gt " << allele_gt.second << endl;
#endif
        if (allele_gt.second > 0) {
            out_variant.alt[allele_gt.second - 1] = allele_gt.first;
        }
        out_variant.alleles[allele_gt.second] = allele_gt.first;

        // update the traversal info
        add_allele_path_to_info(out_variant, allele_gt.second, site_traversals.at(allele_gt.second), false, false); 
    }

    // resolve subpath naming
    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    // in VCF we usually just want a contig
    string contig_name = PathMetadata::parse_locus_name(basepath_name);
    if (contig_name != PathMetadata::NO_LOCUS_NAME) {
        basepath_name = contig_name;
    }
    // fill out the rest of the variant    
    out_variant.sequenceName = basepath_name;
    // +1 to convert to 1-based VCF
    out_variant.position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_offset + 1 + basepath_offset;
    out_variant.id = print_snarl(snarl, false);
    out_variant.filter = "PASS";
    out_variant.updateAlleleIndexes();

    // add the genotype
    out_variant.format.push_back("GT");
    auto& genotype_vector = out_variant.samples[sample_name]["GT"];
    
    stringstream vcf_gt;
    if (!genotype.empty()) {
        for (int i = 0; i < site_genotype.size(); ++i) {
            if (site_genotype[i] == MISSING_ALLELE_MARKER) {
                vcf_gt << ".";
            } else {
                vcf_gt << site_genotype[i];
            }
            if (i != site_genotype.size() - 1) {
                vcf_gt << "/";
            }
        }
    } else {
        for (int i = 0; i < ploidy; ++i) {
            vcf_gt << ".";
            if (i != ploidy - 1) {
                vcf_gt << "/";
            }
        }
    }
                    
    genotype_vector.push_back(vcf_gt.str());

    // add some support info
    snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, call_info, sample_name, out_variant);

    // if genotype_snarls, then we only flatten up to the snarl endpoints
    // (this is when we are in genotyping mode and want consistent calls regardless of the sample)
    int64_t flatten_len_s = 0;
    int64_t flatten_len_e = 0;
    if (genotype_snarls) {
        flatten_len_s = graph.get_length(graph.get_handle(snarl.start().node_id()));
        assert(flatten_len_s >= 0);
        flatten_len_e = graph.get_length(graph.get_handle(snarl.end().node_id()));
    }
    // clean up the alleles to not have so man common prefixes
    flatten_common_allele_ends(out_variant, true, flatten_len_e);
    flatten_common_allele_ends(out_variant, false, flatten_len_s);

    // Merge near-identical called ALT alleles (vg call -L), turning 1/2 into 1/1.  Placed here on
    // purpose: after update_vcf_info so the genotyper saw every candidate, and after flattening so
    // the surviving allele's string, POS and REF are byte-identical to a run without -L (a shorter
    // allele list can share a longer prefix and flatten further).  The missing-allele fixup below
    // still runs after it, and is unaffected: merging is ALT-vs-ALT so it never empties alt.
    merge_similar_alleles(graph, site_traversals, site_genotype, sample_name, out_variant);
#ifdef debug
    for (int i = 0; i < site_traversals.size(); ++i) {
        cerr << " site trav[" << i << "]=" << pb2json(site_traversals[i]) << endl;
    }
    for (int i = 0; i < site_genotype.size(); ++i) {
        cerr << " site geno[" << i << "]=" << site_genotype[i] << endl;
    }
#endif

    // If genotype contains missing allele but no ALT, add * as ALT to emit valid VCF
    // This happens when one parent haplotype doesn't traverse a nested child snarl
    bool has_missing = std::find(site_genotype.begin(), site_genotype.end(), MISSING_ALLELE_MARKER) != site_genotype.end();
    if (has_missing && out_variant.alt.empty()) {
        out_variant.alt.push_back("*");
        out_variant.alleles.push_back("*");
        out_variant.info["AT"].push_back(".");
    }

    if (genotype_snarls || !out_variant.alt.empty()) {
        bool added = add_variant(out_variant);
        if (!added) {
            stringstream ss;
            ss << out_variant;
            cerr << "Warning [vg call]: Skipping variant at " << out_variant.sequenceName << ":" << out_variant.position
                 << " with ID=" << out_variant.id << " because its line length of " << ss.str().length() << " exceeds vg's limit of "
                 << VCFOutputCaller::max_vcf_line_length << endl;
        }
        return added;
    }
    // Same as in Deconstructor::deconstruct_site: a site with nothing to report is still the only
    // thing that knows where its children sit, so keep its reference interval for their RC/RS/RD.
    // Kept in step with deconstruct deliberately -- the two tools should annotate a gref graph the
    // same way -- even though vg call's records are on reference paths, where the old self-reference
    // fallback was at least a real coordinate.
    if (include_nested) {
        suppressed_ref_info[omp_get_thread_num()][out_variant.id] =
            {out_variant.sequenceName, static_cast<size_t>(out_variant.position),
             out_variant.ref.length()};
    }
    return true;
}

tuple<int64_t, int64_t, bool, step_handle_t, step_handle_t> VCFOutputCaller::get_ref_interval(
    const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name) const {
    path_handle_t path_handle = graph.get_path_handle(ref_path_name);

    handle_t start_handle = graph.get_handle(snarl.start().node_id(), snarl.start().backward());
    map<size_t, step_handle_t> start_steps;
    graph.for_each_step_on_handle(start_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                start_steps[graph.get_position_of_step(step)] = step;
            }
        });

    handle_t end_handle = graph.get_handle(snarl.end().node_id(), snarl.end().backward());
    map<size_t, step_handle_t> end_steps;
    graph.for_each_step_on_handle(end_handle, [&](step_handle_t step) {
            if (graph.get_path_handle_of_step(step) == path_handle) {
                end_steps[graph.get_position_of_step(step)] = step;
            }
        });

    assert(start_steps.size() > 0 && end_steps.size() > 0);
    step_handle_t start_step = start_steps.begin()->second;
    step_handle_t end_step = end_steps.begin()->second;
    // just because we found a pair of steps on our path that correspond to the snarl ends, doesn't
    // mean the path threads the snarl.  verify that we can actaully walk, either forwards or backwards
    // along the path from the start node and hit then end node in the right orientation. 
    bool start_rev = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
    bool end_rev = graph.get_is_reverse(graph.get_handle_of_step(end_step)) != snarl.end().backward();
    bool found_end = start_rev == end_rev && start_rev == start_steps.begin()->first > end_steps.begin()->first;
        
    // if we're on a cycle, we keep our start step and find the end step by scanning the path
    if (start_steps.size() > 1 || end_steps.size() > 1) {
        found_end = false;
        // try each start step
        for (auto i = start_steps.begin(); i != start_steps.end() && !found_end; ++i) {
            start_step = i->second;
            bool scan_backward = graph.get_is_reverse(graph.get_handle_of_step(start_step)) != snarl.start().backward();
            if (scan_backward) {
                // if we're going backward, we expect to reach the end backward
                end_handle = graph.get_handle(snarl.end().node_id(), !snarl.end().backward());
            }            
            if (scan_backward) {
                for (step_handle_t cur_step = start_step; graph.has_previous_step(cur_step) && !found_end;
                     cur_step = graph.get_previous_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            } else {
                for (step_handle_t cur_step = start_step; graph.has_next_step(cur_step) && !found_end;
                     cur_step = graph.get_next_step(cur_step)) {
                    if (graph.get_handle_of_step(cur_step) == end_handle) {
                        end_step = cur_step;
                        found_end = true;
                    }
                }
            }
        }
    }
    int64_t start_position = start_steps.begin()->first;
    step_handle_t out_start_step = start_step;
    int64_t end_position = end_step == end_steps.begin()->second ? end_steps.begin()->first : graph.get_position_of_step(end_step);
    step_handle_t out_end_step = end_step == end_steps.begin()->second ? end_steps.begin()->second : end_step;
    bool backward = end_position < start_position;
    

    if (!found_end) {
        // oops, once of the above checks failed.  we tell caller we coudlnt find by hacking in a -1 coordinate.
        start_position = -1;
        end_position = -1;
    }

    if (backward) {
        return make_tuple(end_position, start_position, backward, out_end_step, out_start_step);
    } else {
        return make_tuple(start_position, end_position, backward, out_start_step, out_end_step);
    }
}

pair<string, int64_t> VCFOutputCaller::get_ref_position(const PathPositionHandleGraph& graph, const Snarl& snarl, const string& ref_path_name,
                                                        int64_t ref_path_offset) const {

    subrange_t subrange;
    string basepath_name = Paths::strip_subrange(ref_path_name, &subrange);
    size_t basepath_offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
    // +1 to convert to 1-based VCF
    int64_t position = get<0>(get_ref_interval(graph, snarl, ref_path_name)) + ref_path_offset + 1 + basepath_offset;
    return make_pair(basepath_name, position);
}

void VCFOutputCaller::flatten_common_allele_ends(vcflib::Variant& variant, bool backward, size_t len_override) const {
    if (variant.alt.size() == 0) {
        return;
    }

    // find the minimum allele length to make sure we don't delete an entire allele
    size_t min_allele_len = variant.alleles[0].length();
    for (int i = 1; i < variant.alleles.size(); ++i) {
        min_allele_len = std::min(min_allele_len, variant.alleles[i].length());
    }

    // the maximum number of bases we want ot zip up, applying override if provided
    size_t max_flatten_len = len_override > 0 ? len_override : min_allele_len;
    
    // want to leave at least one in the reference position
    if (max_flatten_len == min_allele_len) {
        --max_flatten_len;
    }
    
    bool match = true;
    int shared_prefix_len = 0;
    for (int i = 0; i < max_flatten_len && match; ++i) {
        char c1 = std::toupper(variant.alleles[0][!backward ? i : variant.alleles[0].length() - 1 - i]);
        for (int j = 1; j < variant.alleles.size() && match; ++j) {
            char c2 = std::toupper(variant.alleles[j][!backward ? i : variant.alleles[j].length() - 1 - i]);
            match = c1 == c2;
        }
        if (match) {
            ++shared_prefix_len;
        }
    }

    if (!backward) {
        variant.position += shared_prefix_len;
    }
    for (int i = 0; i < variant.alleles.size(); ++i) {
        if (!backward) {
            variant.alleles[i] = variant.alleles[i].substr(shared_prefix_len);
        } else {
            variant.alleles[i] = variant.alleles[i].substr(0, variant.alleles[i].length() - shared_prefix_len);
        }
        if (i == 0) {
            variant.ref = variant.alleles[i];
        } else {
            variant.alt[i - 1] = variant.alleles[i];
        }
    }
}

string VCFOutputCaller::nesting_info_headers() {
    stringstream ss;
    ss << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree counting only ancestors whose record is on this record's own reference contig (0=top level for this contig)\">" << endl;
    ss << "##INFO=<ID=CH,Number=1,Type=Integer,Description=\"Number of ancestors in the VCF whose record is on a different reference contig than the one below it, ie nesting steps between VCF reference contigs\">" << endl;
    ss << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    ss << "##INFO=<ID=RC,Number=1,Type=String,Description=\"CHROM of the topmost ancestor record in this VCF, or this record's own CHROM when it has none. On a gref fragment, where that own CHROM would be no use, the enclosing site is named even if it produced no record of its own; the tags are absent when there is no such site either\">" << endl;
    ss << "##INFO=<ID=RS,Number=1,Type=Integer,Description=\"Start of the site named by RC: the POS of its record, or where the site begins when it produced none. A position on that contig, not a span of the snarl, so it can precede the sequence this record describes\">" << endl;
    ss << "##INFO=<ID=RD,Number=1,Type=Integer,Description=\"End of the site named by RC: RS plus the length of that site's REF allele\">" << endl;
    return ss.str();
}

string VCFOutputCaller::print_snarl(const HandleGraph* graph, const handle_t& snarl_start,
                                    const handle_t& snarl_end, bool in_brackets) const {
    Snarl snarl;
    Visit* start = snarl.mutable_start();
    start->set_node_id(graph->get_id(snarl_start));
    start->set_backward(graph->get_is_reverse(snarl_start));
    Visit* end = snarl.mutable_end();
    end->set_node_id(graph->get_id(snarl_end));
    end->set_backward(graph->get_is_reverse(snarl_end));
    return this->print_snarl(snarl, in_brackets);
}

string VCFOutputCaller::print_snarl(const Snarl& snarl, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    nid_t start_node_id = snarl.start().node_id();
    nid_t end_node_id = snarl.end().node_id();
    string start_node = std::to_string(start_node_id);
    string end_node = std::to_string(end_node_id);
    if (translation) {
        auto i = translation->find(start_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + start_node + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node_id);
        if (i == translation->end()) {
            throw runtime_error("Error [VCFOutputCaller]: Unable to find node " + end_node + " in translation file");
        }
        end_node = i->second.first;
    }
    stringstream ss;
    if (in_brackets) {
        ss << "(";
    }
    ss << (snarl.start().backward() ? "<" : ">") << start_node << (snarl.end().backward() ? "<" : ">") << end_node;
    if (in_brackets) {
        ss << ")";
    }
    return ss.str();
}
string VCFOutputCaller::print_flipped_snarl(const Snarl& snarl, bool in_brackets) const {
    // todo, should we canonicalize here by putting lexicographic lowest node first?
    Snarl flipped_snarl;
    flipped_snarl.mutable_start()->set_node_id(snarl.end().node_id());
    flipped_snarl.mutable_start()->set_backward(!snarl.end().backward());
    flipped_snarl.mutable_end()->set_node_id(snarl.start().node_id());
    flipped_snarl.mutable_end()->set_backward(!snarl.start().backward());
    return print_snarl(flipped_snarl, in_brackets);
}

void VCFOutputCaller::scan_snarl(const string& allele_string, function<void(const string&, Snarl&)> callback) const {
    int left = -1;
    int last = 0;
    Snarl snarl;
    string frag;
    for (int i = 0; i < allele_string.length(); ++i) {
        if (allele_string[i] == '(') {
            assert(left == -1);
            if (last < i) {
                frag = allele_string.substr(last, i-last);
                callback(frag, snarl);
            }
            left = i;
        } else if (allele_string[i] == ')') {
            assert(left >= 0 && i > left + 3);
            frag = allele_string.substr(left + 1, i - left - 1);
            auto toks = split_delims(frag, "><");
            assert(toks.size() == 2);
            assert(frag[0] == '<' || frag[0] == '>');
            int64_t start = std::stoi(toks[0]);
            snarl.mutable_start()->set_node_id(start);
            snarl.mutable_start()->set_backward(frag[0] == '<');
            assert(frag[toks[0].size() + 1] == '<' || frag[toks[0].size() + 1] == '>');
            int64_t end = std::stoi(toks[1]);
            snarl.mutable_end()->set_node_id(abs(end));
            snarl.mutable_end()->set_backward(frag[toks[0].size() + 1] == '<');
            callback("", snarl);
            left = -1;
            last = i + 1;
        }
    }
    if (last == 0) {
        callback(allele_string, snarl);
    } else {
        frag = allele_string.substr(last);
        callback(frag, snarl);
    }
}

void VCFOutputCaller::update_nesting_info_tags(const SnarlManager* snarl_manager) {

    // index the snarl tree by name
    unordered_map<string, const Snarl*> name_to_snarl;
    Snarl flipped_snarl;
    snarl_manager->for_each_snarl_preorder([&](const Snarl* snarl) {
            name_to_snarl[print_snarl(*snarl)] = snarl;
            // also add a map from the flipped snarl (as call sometimes messes with orientation)
            flipped_snarl.mutable_start()->set_node_id(snarl->end().node_id());
            flipped_snarl.mutable_start()->set_backward(!snarl->end().backward());
            flipped_snarl.mutable_end()->set_node_id(snarl->start().node_id());
            flipped_snarl.mutable_end()->set_backward(!snarl->start().backward());
            name_to_snarl[print_snarl(flipped_snarl)] = snarl;
        });

    // Merge the per-thread suppressed-site intervals collected during calling.  These are sites
    // that never reached the VCF, so pass 1 below cannot see them, but a record nested under one
    // has no other way to name a reference position.
    unordered_map<string, SuppressedRef> suppressed_ref;
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
    for (auto& thread_buf : output_variants) {
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 4);
            chrom_of_name.emplace(toks[2], intern_chrom(toks[0]));
        }
    }

    // pass 2) identify top-level snarls (those with no ancestors in VCF)
    // and store reference info only for them
    struct RefInfo {
        string chrom;
        size_t pos;
        size_t ref_len;
    };
    // Keyed by snarl name, then by (chrom, pos), because a snarl ID can carry more than one
    // record: a cyclic reference emits two, both with the same ID (see
    // nesting/cyclic_ref_multiple_variants.gfa, which gives two <5<1 records at POS 20 and 44).
    // A plain name -> RefInfo map was last-write-wins, so both records were handed the
    // surviving one's interval and the record at POS 20 reported RS=44.  The inner map is
    // ordered so that picking begin() is deterministic regardless of thread scheduling.
    unordered_map<string, map<pair<string, size_t>, size_t>> top_level_ref_info;

    // Helper to check if a snarl is top-level (no ancestors in VCF)
    auto is_top_level = [&](const string& name) -> bool {
        auto it = name_to_snarl.find(name);
        if (it == name_to_snarl.end()) return true; // not found, treat as top-level
        const Snarl* snarl = it->second;
        while ((snarl = snarl_manager->parent_of(snarl))) {
            string cur_name = print_snarl(*snarl);
            string flipped_name = print_flipped_snarl(*snarl);
            if (chrom_of_name.count(cur_name) || chrom_of_name.count(flipped_name)) {
                return false; // has ancestor in VCF
            }
        }
        return true; // no ancestors in VCF
    };

    // Second pass through variants to extract ref info only for top-level snarls
    for (auto& thread_buf : output_variants) {
        for (auto& output_variant_record : thread_buf) {
            string output_variant_string;
            int ret = zstdutil::DecompressString(output_variant_record.second, output_variant_string);
            assert(ret == 0);
            vector<string> toks = split_delims(output_variant_string, "\t", 5);
            const string& name = toks[2];
            if (is_top_level(name)) {
                top_level_ref_info[name][make_pair(toks[0], static_cast<size_t>(stoul(toks[1])))] =
                    toks[3].length();
            }
        }
    }

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
        const Snarl* snarl = name_to_snarl.at(name);

        assert(snarl != nullptr);
        // walk up the snarl tree
        while ((snarl = snarl_manager->parent_of(snarl))) {
            string cur_name = print_snarl(*snarl);

            // Since it is possible that the snarl is actually flipped in the vcf, check for the flipped version too
            string flipped_name = print_flipped_snarl(*snarl);
            const string* hit = nullptr;
            if (chrom_of_name.count(cur_name)) {
                // only count snarls that are in the vcf
                hit = &cur_name;
            } else if (chrom_of_name.count(flipped_name)) {
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
            //string& output_variant_string = output_variant_record.second;
            vector<string> toks = split_delims(output_variant_string, "\t", 9);
            const string& name = toks[2];

            auto [contig_level, contig_hops, parent_name, top_level_name, ref_chrom_name,
                  suppressed_name] = get_nesting_tags(name, toks[0]);
            // LV is the level within this record's own reference contig.  It used to count
            // ancestors across every contig: for a VCF with a single reference contig the two
            // are identical, but once gref fragments give the insides of insertions their own
            // contigs, the whole-file count is not what a level filter wants.  No gref-contig
            // record was ever at LV=0 under the old definition, so `vcfbub -l 0` deleted every
            // one of them.
            string nesting_tags = ";LV=" + std::to_string(contig_level);
            nesting_tags += ";CH=" + std::to_string(contig_hops);
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
            RefInfo top_ref;
            if (ref_source == nullptr) {
                if (!GrefCover::is_gref_name(toks[0])) {
                    // Not on a gref fragment, so our own interval is already a position a reader
                    // can look up, and it is the narrower answer of the two.  Keep it rather than
                    // reach for an enclosing site that produced no record: doing that would
                    // repoint every such record on a reference contig at a site LV and CH say it
                    // has no ancestor in.  The records that need the reach are the fragments
                    // below, which have no usable coordinate of their own.
                    top_ref = {toks[0], static_cast<size_t>(stoul(toks[1])), toks[3].length()};
                } else {
                    auto sup_it = suppressed_name.empty() ? suppressed_ref.end()
                                                          : suppressed_ref.find(suppressed_name);
                    if (sup_it != suppressed_ref.end()) {
                        top_ref = {sup_it->second.chrom, sup_it->second.pos, sup_it->second.ref_len};
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
                nesting_tags += ";RC=" + top_ref.chrom;
                nesting_tags += ";RS=" + std::to_string(top_ref.pos);
                nesting_tags += ";RD=" + std::to_string(top_ref.pos + top_ref.ref_len);
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

