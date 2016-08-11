#include "readfilter.hpp"
#include "IntervalTree.h"

#include <fstream>
#include <sstream>

namespace vg {

using namespace std;

bool ReadFilter::has_repeat(Alignment& aln, int k) {
    if (k == 0) {
        return false;
    }
    const string& s = aln.sequence();
    for (int i = 1; i <= 2 * k; ++i) {
        int covered = 0;
        bool ffound = true;
        bool bfound = true;
        for (int j = 1; (ffound || bfound) && (j + 1) * i < s.length(); ++j) {
            ffound = ffound && s.substr(0, i) == s.substr(j * i, i);
            bfound = bfound && s.substr(s.length() - i, i) == s.substr(s.length() - i - j * i, i);
            if (ffound || bfound) {
                covered += i;
            }
        }
        if (covered >= k) {
            return true;
        }
    }
    return false;
}

int ReadFilter::filter(istream* alignment_stream) {

    // name helper for output
    function<string(int)> chunk_name = [this](int num) -> string {
        stringstream ss;
        ss << outbase << "-" << num << ".gam";
        return ss.str();
    };

    // index regions by their inclusive ranges
    vector<Interval<int, int64_t> > interval_list;
    vector<Region> regions;
    // use strings instead of ofstreams because worried about too many handles
    vector<string> chunk_names;
    vector<bool> chunk_new; // flag if write or append

    // parse a bed, for now this is only way to do regions.  note
    // this operation converts from 0-based BED to 1-based inclusive VCF
    if (!regions_file.empty()) {
        if (outbase.empty()) {
            cerr << "-B option required with -R" << endl;
            exit(1);
        }
        parse_bed_regions(regions_file, regions);
        if (regions.empty()) {
            cerr << "No regions read from BED file, doing whole graph" << endl;
        }
    }

    // empty region, do everything
    if (regions.empty()) {
        // we handle empty intervals as special case when looking up, otherwise,
        // need to insert giant interval here.
        chunk_names.push_back(outbase.empty() ? "-" : chunk_name(0));
    }

    // otherwise, need to extract regions with xg
    else {
        if (xg_name.empty()) {
            cerr << "xg index required for -R option" << endl;
            exit(1);
        }
        // read the xg index
        ifstream xg_stream(xg_name);
        if(!xg_stream) {
            cerr << "Unable to open xg index: " << xg_name << endl;
            exit(1);
        }
        xg::XG xindex(xg_stream);

        // fill in the map using xg index
        // relies entirely on the assumption that are path chunks
        // are perfectly spanned by an id range
        for (int i = 0; i < regions.size(); ++i) {
            Graph graph;
            int rank = xindex.path_rank(regions[i].seq);
            int path_size = rank == 0 ? 0 : xindex.path_length(regions[i].seq);

            if (regions[i].start >= path_size) {
                cerr << "Unable to find region in index: " << regions[i].seq << ":" << regions[i].start
                     << "-" << regions[i].end << endl;
            } else {
                // clip region on end of path
                regions[i].end = min(path_size - 1, regions[i].end);
                // do path node query
                // convert to 0-based coordinates as this seems to be what xg wants
                xindex.get_path_range(regions[i].seq, regions[i].start - 1, regions[i].end - 1, graph);
                if (context_size > 0) {
                    xindex.expand_context(graph, context_size);
                }
            }
            // find node range of graph, without bothering to build vg indices..
            int64_t min_id = numeric_limits<int64_t>::max();
            int64_t max_id = 0;
            for (int j = 0; j < graph.node_size(); ++j) {
                min_id = min(min_id, (int64_t)graph.node(j).id());
                max_id = max(max_id, (int64_t)graph.node(j).id());
            }
            // map the chunk id to a name
            chunk_names.push_back(chunk_name(i));

            // map the node range to the chunk id.
            if (graph.node_size() > 0) {
                interval_list.push_back(Interval<int, int64_t>(min_id, max_id, i));
                assert(chunk_names.size() == i + 1);
            }
        }
    }

    // index chunk regions
    IntervalTree<int, int64_t> region_map(interval_list);

    // which chunk(s) does a gam belong to?
    function<void(Alignment&, vector<int>&)> get_chunks = [&region_map, &regions](Alignment& aln, vector<int>& chunks) {
        // speed up case where no chunking
        if (regions.empty()) {
            chunks.push_back(0);
        } else {
            int64_t min_aln_id = numeric_limits<int64_t>::max();
            int64_t max_aln_id = -1;
            for (int i = 0; i < aln.path().mapping_size(); ++i) {
                const Mapping& mapping = aln.path().mapping(i);
                min_aln_id = min(min_aln_id, (int64_t)mapping.position().node_id());
                max_aln_id = max(max_aln_id, (int64_t)mapping.position().node_id());
            }
            vector<Interval<int, int64_t> > found_ranges;
            region_map.findOverlapping(min_aln_id, max_aln_id, found_ranges);
            for (auto& interval : found_ranges) {
                chunks.push_back(interval.value);
            }
        }
    };

    // buffered output (one buffer per chunk)
    vector<vector<Alignment> > buffer(chunk_names.size());
    int cur_buffer = -1;
    static const int buffer_size = 1000; // we let this be off by 1
    function<Alignment&(uint64_t)> write_buffer = [&buffer, &cur_buffer](uint64_t i) -> Alignment& {
        return buffer[cur_buffer][i];
    };
    // remember if write or append
    vector<bool> chunk_append(chunk_names.size(), false);

    // flush a buffer specified by cur_buffer to target in chunk_names, and clear it
    function<void()> flush_buffer = [&buffer, &cur_buffer, &write_buffer, &chunk_names, &chunk_append]() {
        ofstream outfile;
        auto& outbuf = chunk_names[cur_buffer] == "-" ? cout : outfile;
        if (chunk_names[cur_buffer] != "-") {
            outfile.open(chunk_names[cur_buffer], chunk_append[cur_buffer] ? ios::app : ios_base::out);
            chunk_append[cur_buffer] = true;
        }
        stream::write(outbuf, buffer[cur_buffer].size(), write_buffer);
        buffer[cur_buffer].clear();
    };

    // add alignment to all appropriate buffers, flushing as necessary
    // (note cur_buffer variable used here as a hack to specify which buffer is written to)
    function<void(Alignment&)> update_buffers = [&buffer, &cur_buffer, &region_map,
                                                 &write_buffer, &get_chunks, &flush_buffer](Alignment& aln) {
        vector<int> aln_chunks;
        get_chunks(aln, aln_chunks);
        for (auto chunk : aln_chunks) {
            buffer[chunk].push_back(aln);
            if (buffer[chunk].size() >= buffer_size) {
                // flush buffer
                cur_buffer = chunk;
                flush_buffer();
            }
        }
    };

    // keep track of how many reads were dropped by which option
    size_t pri_read_count = 0;
    size_t sec_read_count = 0;
    size_t sec_filtered_count = 0;
    size_t pri_filtered_count = 0;
    size_t min_sec_count = 0;
    size_t min_pri_count = 0;
    size_t min_sec_delta_count = 0;
    size_t min_pri_delta_count = 0;
    size_t max_sec_overhang_count = 0;
    size_t max_pri_overhang_count = 0;
    size_t min_sec_mapq_count = 0;
    size_t min_pri_mapq_count = 0;
    size_t repeat_sec_count = 0;
    size_t repeat_pri_count = 0;

    // for deltas, we keep track of last primary
    Alignment prev;
    bool prev_primary = false;
    bool keep_prev = true;
    double prev_score;

    // quick and dirty filter to see if removing reads that can slip around
    // and still map perfectly helps vg call.  returns true if at either
    // end of read sequence, at least k bases are repetitive, checking repeats
    // of up to size 2k
    
        
    // we assume that every primary alignment has 0 or 1 secondary alignment
    // immediately following in the stream
    function<void(Alignment&)> lambda = [&](Alignment& aln) {
        bool keep = true;
        double score = (double)aln.score();
        double denom = 2. * aln.sequence().length();
        // toggle substitution score
        if (sub_score == true) {
            // hack in ident to replace old counting logic.
            score = aln.identity() * aln.sequence().length();
            denom = aln.sequence().length();
            assert(score <= denom);
        }
        // toggle absolute or fractional score
        if (frac_score == true) {
            if (denom > 0.) {
                score /= denom;
            }
            else {
                assert(score == 0.);
            }
        }
        // compute overhang
        int overhang = 0;
        if (aln.path().mapping_size() > 0) {
            const auto& left_mapping = aln.path().mapping(0);
            if (left_mapping.edit_size() > 0) {
                overhang = left_mapping.edit(0).to_length() - left_mapping.edit(0).from_length();
            }
            const auto& right_mapping = aln.path().mapping(aln.path().mapping_size() - 1);
            if (right_mapping.edit_size() > 0) {
                const auto& edit = right_mapping.edit(right_mapping.edit_size() - 1);
                overhang = max(overhang, edit.to_length() - edit.from_length());
            }
        } else {
            overhang = aln.sequence().length();
        }

        if (aln.is_secondary()) {
            ++sec_read_count;
            assert(prev_primary && aln.name() == prev.name());
            double delta = prev_score - score;
            if (frac_delta == true) {
                delta = prev_score > 0 ? score / prev_score : 0.;
            }

            // filter (current) secondary
            keep = true;
            if (score < min_secondary) {
                ++min_sec_count;
                keep = false;
            }
            if ((keep || verbose) && delta < min_sec_delta) {
                ++min_sec_delta_count;
                keep = false;
            }
            if ((keep || verbose) && overhang > max_overhang) {
                ++max_sec_overhang_count;
                keep = false;
            }
            if ((keep || verbose) && aln.mapping_quality() < min_mapq) {
                ++min_sec_mapq_count;
                keep = false;
            }
            if ((keep || verbose) && has_repeat(aln, repeat_size)) {
                ++repeat_sec_count;
                keep = false;
            }
            if (!keep) {
                ++sec_filtered_count;
            }

            // filter unfiltered previous primary
            if (keep_prev && delta < min_pri_delta) {
                ++min_pri_delta_count;
                ++pri_filtered_count;
                keep_prev = false;
            }
            // add to write buffer
            if (keep) {
                update_buffers(aln);
            }
            if (keep_prev) {
                update_buffers(prev);
            }

            // forget last primary
            prev_primary = false;
            prev_score = -1;
            keep_prev = false;

        } else {
            // this awkward logic where we keep the primary and write in the secondary
            // is because we can only stream one at a time with for_each, but we need
            // to look at pairs (if secondary exists)...
            ++pri_read_count;
            if (keep_prev) {
                update_buffers(prev);
            }

            prev_primary = true;
            prev_score = score;
            prev = aln;
            keep_prev = true;
            if (score < min_primary) {
                ++min_pri_count;
                keep_prev = false;
            }
            if ((keep_prev || verbose) && overhang > max_overhang) {
                ++max_pri_overhang_count;
                keep_prev = false;
            }
            if ((keep_prev || verbose) && aln.mapping_quality() < min_mapq) {
                ++min_pri_mapq_count;
                keep_prev = false;
            }
            if ((keep_prev || verbose) && has_repeat(aln, repeat_size)) {
                ++repeat_pri_count;
                keep_prev = false;
            }
            if (!keep_prev) {
                ++pri_filtered_count;
            }
        }
    };
    stream::for_each(*alignment_stream, lambda);

    // flush buffer if trailing primary to keep
    if (keep_prev) {
        update_buffers(prev);
    }

    for (int i = 0; i < buffer.size(); ++i) {
        if (buffer[i].size() > 0) {
            cur_buffer = i;
            flush_buffer();
        }
    }

    if (verbose) {
        size_t tot_reads = pri_read_count + sec_read_count;
        size_t tot_filtered = pri_filtered_count + sec_filtered_count;
        cerr << "Total Filtered (primary):          " << pri_filtered_count << " / "
             << pri_read_count << endl
             << "Total Filtered (secondary):        " << sec_filtered_count << " / "
             << sec_read_count << endl
             << "Min Identity Filter (primary):     " << min_pri_count << endl
             << "Min Identity Filter (secondary):   " << min_sec_count << endl
             << "Min Delta Filter (primary):        " << min_pri_delta_count << endl
             << "Min Delta Filter (secondary):      " << min_sec_delta_count << endl
             << "Max Overhang Filter (primary):     " << max_pri_overhang_count << endl
             << "Max Overhang Filter (secondary):   " << max_sec_overhang_count << endl
             << "Repeat Ends Filter (primary):     " << repeat_pri_count << endl
             << "Repeat Ends Filter (secondary):   " << repeat_sec_count << endl

             << endl;
    }

    return 0;

}

}
