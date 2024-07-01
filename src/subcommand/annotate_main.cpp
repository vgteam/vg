#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../mapper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../alignment.hpp"
#include "../annotation.hpp"
#include "../gff_reader.hpp"
#include "../region_expander.hpp"
#include "../algorithms/alignment_path_offsets.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

#include "progress_bar.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_annotate(char** argv) {
    cerr << "usage: " << argv[0] << " annotate [options] >output.{gam,vg,tsv}" << endl
         << "graph annotation options:" << endl
         << "    -x, --xg-name FILE     xg index or graph to annotate (required)" << endl
         << "    -b, --bed-name FILE    a BED file to convert to GAM. May repeat." << endl
         << "    -f, --gff-name FILE    a GFF3 file to convert to GAM. May repeat." << endl
         << "    -g, --ggff             output at GGFF subgraph annotation file instead of GAM (requires -s)" << endl
         << "    -s, --snarls FILE      file containing snarls to expand GFF intervals into" << endl
         << "alignment annotation options:" << endl
         << "    -a, --gam FILE         file of Alignments to annotate (required)" << endl
         << "    -x, --xg-name FILE     xg index of the graph against which the Alignments are aligned (required)" << endl
         << "    -p, --positions        annotate alignments with reference positions" << endl
         << "    -m, --multi-position   annotate alignments with multiple reference positions" << endl
         << "    -l, --search-limit N   when annotating with positions, search this far for paths, or -1 to not search (default: 0 (auto from read length))" << endl
         << "    -b, --bed-name FILE    annotate alignments with overlapping region names from this BED. May repeat." << endl
         << "    -n, --novelty          output TSV table with header describing how much of each Alignment is novel" << endl
         << "    -P, --progress         show progress" << endl
         << "    -t, --threads          use the specified number of threads" << endl;
}

/// Find the region of the Mapping's node used by the Mapping, in forward strand space, as start to past_end.
static pair<size_t, size_t> mapping_to_range(const HandleGraph* xg_index, const Mapping& mapping) {
    // How much of the node does it cover?
    auto mapping_length = mapping_from_length(mapping);
    
    // Work out where the start and past-end positions on the node's forward strand are.
    pair<size_t, size_t> node_range;
    if (mapping.position().is_reverse()) {
        // On the reverse strand we need the node length
        // TODO: getting it can be slow
        auto node_length = xg_index->get_length(xg_index->get_handle(mapping.position().node_id()));
        
        node_range.first = node_length - mapping.position().offset() - mapping_length;
        node_range.second = node_length - mapping.position().offset();
    } else {
        // On the forward strand this is easy
        node_range.first = mapping.position().offset();
        node_range.second = node_range.first + mapping_length;
    }
    
    return node_range;
}

using feature_t = pair<pair<size_t, size_t>, const string*>;

/// Look up all ranges in the vector of features that overlap the start to past-end search interval.
/// Return their unique string* interned feature names.
static unordered_set<const string*> find_overlapping(const vector<feature_t>& ranges, pair<size_t, size_t> search_interval) {
    // TODO: We expect at most a few features per node, so just do a linear scan
    
    unordered_set<const string*> to_return;
    
    for (auto& feature : ranges) {
        auto& feature_range = feature.first;
        
        if (feature_range.second > search_interval.first && feature_range.first < search_interval.second) {
            // There is overlap
            to_return.insert(feature.second);
        }
    }
    
    return to_return;
}

int main_annotate(int argc, char** argv) {
    
    if (argc == 2) {
        help_annotate(argv);
        return 1;
    }

    string xg_name;
    vector<string> bed_names;
    vector<string> gff_names;
    string gam_name;
    bool add_positions = false;
    bool add_multiple_positions = false;
    int64_t search_limit = 0;
    bool novelty = false;
    bool output_ggff = false;
    string snarls_name;
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"gam", required_argument, 0, 'a'},
            {"positions", no_argument, 0, 'p'},
            {"multi-positions", no_argument, 0, 'm'},
            {"search-limit", required_argument, 0, 'l'},
            {"xg-name", required_argument, 0, 'x'},
            {"bed-name", required_argument, 0, 'b'},
            {"gff-name", required_argument, 0, 'f'},
            {"ggff", no_argument, 0, 'g'},
            {"snarls", required_argument, 0, 's'},
            {"novelty", no_argument, 0, 'n'},
            {"progress", no_argument, 0, 'P'},
            {"threads", required_argument, 0, 't'},
            {"help", required_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:a:pml:b:f:gs:nt:Ph",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'x':
            xg_name = optarg;
            break;

        case 'a':
            gam_name = optarg;
            break;

        case 'b':
            bed_names.push_back(optarg);
            break;

        case 'f':
            gff_names.push_back(optarg);
            break;
            
        case 'g':
            output_ggff = true;
            break;
                
        case 's':
            snarls_name = optarg;
            break;

        case 'p':
            add_positions = true;
            break;
            
        case 'm':
            add_positions = true;
            add_multiple_positions = true;
            break;
            
        case 'l':
            search_limit = parse<int64_t>(optarg);
            break;
            
        case 'n':
            novelty = true;
            break;
            
        case 't':
            omp_set_num_threads(parse<size_t>(optarg));
            break;

        case 'P':
            show_progress = true;
            break;

        case 'h':
        case '?':
            help_annotate(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }
    
    PathPositionHandleGraph* xg_index = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::ReferencePathOverlayHelper overlay_helper;

    if (!xg_name.empty()) {
        // Read in the XG index
        if (show_progress) {
            std::cerr << "Load graph" << std::endl;
        }
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        if (show_progress) {
            std::cerr << "Apply overlay" << std::endl;
        }
        xg_index = overlay_helper.apply(path_handle_graph.get());
    } else {
        cerr << "error [vg annotate]: no xg index provided" << endl;
        return 1;
    }
    
    
    unique_ptr<SnarlManager> snarl_manager = nullptr;
    if (!snarls_name.empty()) {
        if (show_progress) {
            std::cerr << "Load snarls" << std::endl;
        }
        get_input_file(snarls_name, [&](istream& snarl_stream) {
            snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_stream);
        });
    }
    
    Mapper mapper(xg_index, nullptr, nullptr);
    
    if (!gam_name.empty()) {
        vector<Alignment> buffer;
        
        if (novelty) {
            // We are computing a novelty table.
            // TODO: refactor this into novelty annotation and annotation-to-table conversion.
            if (add_positions || !bed_names.empty()) {
                // We can't make the TSV and also annotate the reads
                cerr << "error [vg annotate]: Cannot annotate reads while computing novelty table" << endl;
                return 1;
            }
            
            cout << "name\tlength.bp\tunaligned.bp\tknown.nodes\tknown.bp\tnovel.nodes\tnovel.bp" << endl;
            function<void(Alignment&)> lambda = [&](Alignment& aln) {
                // count the number of positions in the alignment that aren't in the graph
                int total_bp = aln.sequence().size();
                int unaligned_bp = 0;
                int known_nodes = 0;
                int known_bp = 0;
                int novel_nodes = 0;
                int novel_bp = 0;
                for (auto& mapping : aln.path().mapping()) {
                    if (mapping.has_position()) {
                        auto& pos = mapping.position();
                        if (xg_index->has_node(pos.node_id())) {
                            ++known_nodes;
                            known_bp += mapping_to_length(mapping);
                        } else {
                            ++novel_nodes;
                            novel_bp += mapping_to_length(mapping);
                        }
                    } else {
                        unaligned_bp += mapping_to_length(mapping);
                    }
                }
                cout << aln.name() << "\t"
                << total_bp << "\t"
                << unaligned_bp << "\t"
                << known_nodes << "\t"
                << known_bp << "\t"
                << novel_nodes << "\t"
                << novel_bp << endl;
            };
            get_input_file(gam_name, [&](istream& in) {
                vg::Progressive::with_progress(show_progress, "Read reads", [&](const std::function<void(size_t, size_t)>& progress) {
                    vg::io::for_each(in, lambda, progress);
                });
            });
        } else {
            // We are annotating the actual reads
            
            // Make per-thread buffers for writing them
            vector<vector<Alignment>> buffers;
            buffers.resize(vg::get_thread_count());
            
            // We will need to track mappings from graph node regions to BED features.
            // We don't want each of those mappings to have a copy of the feature name, because that could be big.
            // So we intern the feature name strings.
            
            // This will map each string to its cannonical copy.
            // TODO: this ends up storing two copies, but is easier than writing the map logic on the unique_ptr value.
            unordered_map<string, unique_ptr<string>> feature_names;
            auto intern = [&feature_names](const string& value) -> const string* {
                if (!feature_names.count(value)) {
                    // If it isn't in the map, put it in
                    feature_names[value] = unique_ptr<string>(new string(value));
                    
                }
                return feature_names[value].get();
            };
            
            // This will hold, for each graph node, the start to past-end
            // regions occupied by BED features, and the names of those
            // features.
            // TODO: To enable better than linear search within a node, we ought to sort these.
            unordered_map<vg::id_t, vector<feature_t>> features_on_node;
            // TODO: We can't use a Paths object because we lack unique names 
            
            for (auto& bed_name : bed_names) {
                // If there are BED files, load them up
                
                get_input_file(bed_name, [&](istream& bed_stream) {
                    // Load all the BED regions as Alignments embedded in the graph.
                    vector<Alignment> bed_regions;
                    parse_bed_regions(bed_stream, xg_index, &bed_regions);
                    
                    for (auto& region : bed_regions) {
                        // For each region in the BED
                        
                        // Get the cannonical copy of its name (which may be "")
                        const string* interned_name = intern(region.name());
                        
                        for (auto& mapping : region.path().mapping()) {
                            // Scan the Mappings. We know each Mapping will be all perfect matches.
                            
                            // Record that the alignment covers the given region on the given node.
                            features_on_node[mapping.position().node_id()].emplace_back(mapping_to_range(xg_index, mapping),
                                interned_name);
                        }
                    }
                });
            }
            
            get_input_file(gam_name, [&](istream& in) {
                vg::Progressive::with_progress(show_progress, "Read reads", [&](const std::function<void(size_t, size_t)>& progress) {
                    vg::io::for_each_parallel<Alignment>(in, [&](Alignment& aln) {
                        // For each read
                        
                        if (add_positions) {
                            // Annotate it with its initial position on each path it touches
                            aln.clear_refpos();
                            if (add_multiple_positions) {
                                // One position per node
                                vg::algorithms::annotate_with_node_path_positions(*mapper.xindex, aln, search_limit);
                            } else {
                                // One position per alignment
                                vg::algorithms::annotate_with_initial_path_positions(*mapper.xindex, aln, search_limit);
                            }
                        }
                        
                        if (!features_on_node.empty()) {
                            // We want to annotate with BED feature overlaps as well.
                            unordered_set<const string*> touched_features;
                            
                            for (auto& mapping : aln.path().mapping()) {
                                // For each mapping
                                
                                auto node_id = mapping.position().node_id();
                                auto features = features_on_node.find(node_id);
                                if (features != features_on_node.end()) {
                                    // Some things occur on this node. Find the overlaps with the part of the node touched by this read.
                                    auto overlapping = find_overlapping(features->second, mapping_to_range(xg_index, mapping));
                                    // Save them all to the set (to remove duplicates)
                                    copy(overlapping.begin(), overlapping.end(), inserter(touched_features, touched_features.begin()));
                                }
                            }
                            
                            // Convert the string pointers to actual string copies, for annotation API.
                            // Make sure to use an ordered set here to sort, to make output deterministic.
                            set<string> feature_names;
                            for (const string* name : touched_features) {
                                feature_names.insert(*name);
                            }
                            
                            // Annotate the read with the feature name strings.
                            set_annotation(aln, "features", feature_names);
                        }
                        
                        // Output the alignment
                        auto& buffer = buffers.at(omp_get_thread_num());
                        buffer.emplace_back(std::move(aln));
                        vg::io::write_buffered(cout, buffer, 1000);
                    }, 256, progress);
                });
            });
        
            for (auto& buffer : buffers) {
                // Finish each buffer
                vg::io::write_buffered(cout, buffer, 0);
            }
        }
    }
    else {
        // Annotating the graph. We must do something.
        if (bed_names.empty() && gff_names.empty()) {
            // We weren't asked to do anything.
            cerr << "error [vg annotate]: only GAM, BED, or GFF3/GTF annotation is implemented" << endl;
            return 1;
        }
    
        if (output_ggff) {
            
            if (!bed_names.empty()) {
                cerr << "error [vg annotate] BED conversion to GGFF is not currently supported. Convert to GFF3 first." << endl;
                return 1;
            }
            
            // define a function that converts to GGFF
            RegionExpander region_expander(&(*xg_index), &(*snarl_manager));
            function<void(const GFFRecord&)> output_ggff_record = [&](const GFFRecord& record) {
                
                auto subgraph = region_expander.expanded_subgraph(record);
                
                if (subgraph.empty()) {
                    cout << ".";
                }
                
                for (auto iter = subgraph.begin(); iter != subgraph.end(); ) {
                    
                    cout << iter->first.first << "[" << iter->second.first << ":" << iter->second.second << "]";
                    if (iter->first.second) {
                        cout << "-";
                    }
                    else {
                        cout << "+";
                    }
                    
                    ++iter;
                    if (iter != subgraph.end()) {
                        cout << ",";
                    }
                }
                cout << "\t";
                
                if (record.source.empty()) {
                    cout << ".";
                }
                else {
                    cout << record.source;
                }
                cout << "\t";
                
                if (record.type.empty()) {
                    cout << ".";
                }
                else {
                    cout << record.type;
                }
                cout << "\t";
                
                if (isnan(record.score)) {
                    cout << ".";
                }
                else {
                    cout << record.score;
                }
                cout << "\t";
                
                if (record.phase == -1) {
                    cout << ".";
                }
                else {
                    cout << record.phase;
                }
                cout << "\t";
                
                if (record.attributes.empty()) {
                    cout << ".";
                }
                else {
                    cout << record.attributes;
                }
                cout << "\n";
            };
            
            for (auto& gff_name : gff_names) {
                get_input_file(gff_name, [&](istream& gff_stream) {
                    GFFReader gff_reader(gff_stream);
                    gff_reader.for_each_gff_record(output_ggff_record);
                });
            }
        }
        else {
            for (auto& bed_name : bed_names) {
                // Convert each BED file to GAM
                get_input_file(bed_name, [&](istream& bed_stream) {
                    vector<Alignment> buffer;
                    parse_bed_regions(bed_stream, xg_index, &buffer);
                    vg::io::write_buffered(cout, buffer, 0); // flush
                });
                
                // TODO: We'll get an EOF marker per input file.
            }
            
            for (auto& gff_name : gff_names) {
                get_input_file(gff_name, [&](istream& gff_stream) {
                    vector<Alignment> buffer;
                    parse_gff_regions(gff_stream, xg_index, &buffer);
                    vg::io::write_buffered(cout, buffer, 0); // flush
                });
            }
        }
    }

    return 0;
}

static Subcommand vg_annotate("annotate", "annotate alignments with graphs and graphs with alignments",
                              main_annotate);
