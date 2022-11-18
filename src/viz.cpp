#include "viz.hpp"
#include <regex>

namespace vg {

Viz::Viz(PathHandleGraph* x, vector<Packer>* p, const vector<string>& n, const string& o, int w, int h, bool c, bool d, bool t) {
    init(x, p, n, o, w, h, c, d, t);
}

void Viz::init(PathHandleGraph* x, vector<Packer>* p, const vector<string>& n, const string& o, int w, int h, bool c, bool d, bool t) {
    xgidx = x;
    packs = p;
    pack_names = n;
    outfile = o;
    show_cnv = c;
    show_dna = d;
    show_paths = t;
    uint64_t rank = 0;
    xgidx->for_each_handle([&](const handle_t& handle) {
            id_rank_map[xgidx->get_id(handle)] = rank++;
            seq_length += xgidx->get_length(handle);
        });
    compute_borders_and_dimensions();
    /*
    left_border = 8;
    top_border = 8;
    image_height = (h ? h : rendered_height() + top_border*2);
    image_width = (w ? w : xgidx->seq_length + xgidx->node_count + left_border*2);
    */
    std::regex svgbase(".svg$");
    std::regex pngbase(".png$");
    if (std::regex_search(outfile, svgbase)) {
        output_svg = true;
        surface = cairo_svg_surface_create(outfile.c_str(), image_width, image_height);
    } else {
        surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, image_width, image_height);
    }
    if (std::regex_search(outfile, pngbase)) {
        output_png = true;
    }
    cr = cairo_create(surface);
}

void Viz::compute_borders_and_dimensions(void) {
    top_border = 4;
    xgidx->for_each_handle([&](const handle_t& h) {
            id_t id = xgidx->get_id(h);
            double s = node_offset(id);
            size_t l = xgidx->get_length(xgidx->get_handle(id));
            xgidx->follow_edges(h, false, [&](const handle_t& o) {
                    id_t id2 = xgidx->get_id(o);
                    double s2 = node_offset(id2);
                    double x = s+l;
                    int delta = s2 - x;
                    double w = pow(log(abs(delta)+1), 1.5);
                    int xdiff = (delta < 0 ? -w : w)/2;
                    int ydiff = w*2;
                    top_border = max(ydiff, top_border);
                    return true;
                });
        });
    int height = top_border + 4;
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 0, 0);
    cr = cairo_create(surface);
    if (show_paths) {
        xgidx->for_each_path_handle([&](const path_handle_t& path) {
                string path_name = xgidx->get_path_name(path);
                cairo_text_extents_t te;
                cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
                cairo_set_font_size(cr, 1);
                cairo_text_extents(cr, path_name.c_str(), &te);
                left_border = max(left_border, (int)round(te.width+2));
                height += 2;
            });
    }
    for (int i = 0; i < packs->size(); ++i) {
        auto& pack_name = pack_names[i];
        cairo_text_extents_t te;
        cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_set_font_size(cr, 1);
        cairo_text_extents(cr, pack_name.c_str(), &te);
        left_border = max(left_border, (int)round(te.width+2));
        height += 2;
    }
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    cr = nullptr;
    surface = nullptr;
    height += 2;
    image_width = seq_length + xgidx->get_node_count() + left_border * 2;
    image_height = height;
}

double Viz::node_offset(id_t id) {
    // provides 1 extra pixel of space per node plus border
    return dynamic_cast<VectorizableHandleGraph*>(xgidx)->node_vector_offset(id) + id_to_rank(id) + left_border;
}

double Viz::nodes_before_offset(size_t pos) {
    // how many nodes have there been before the current position in the seq vector
    return id_to_rank(dynamic_cast<VectorizableHandleGraph*>(xgidx)->node_at_vector_offset(pos+1));
}

uint64_t Viz::id_to_rank(nid_t id) {
    return id_rank_map[id];
}

tuple<double, double, double> hash_to_rgb(const string& str, double min_sum) {
    std::hash<std::string> hash_fn;
    std::size_t h = hash_fn(str);
    uint16_t x, y, z;
    // use 48 bits
    x = h >> 16;
    y = h >> 32;
    z = h >> 48;
    double r, g, b;
    r = (double)x/(double)numeric_limits<uint16_t>::max();
    g = (double)y/(double)numeric_limits<uint16_t>::max();
    b = (double)z/(double)numeric_limits<uint16_t>::max();
    if (r+g+b < min_sum) {
        double to_add = min_sum/3;
        r += to_add;
        g += to_add;
        b += to_add;
    }
    //cerr << "hash " << str << " " << h << " --> " << x << " " << y << " " << z << " --> " << r << " " << g << " " << b << endl;
    return make_tuple(r, g, b);
}

void Viz::set_hash_color(const string& str) {
    auto c = hash_to_rgb(str, 0.5);
    cairo_set_source_rgb(cr, get<0>(c), get<1>(c), get<2>(c));
}

void Viz::draw_graph(void) {
    cairo_set_source_rgb(cr, 0, 0, 0);
    int y_pos = top_border;
    xgidx->for_each_handle([&](const handle_t& h) {
            // get the start and end position of the node relative to the image
            id_t id = xgidx->get_id(h);
            double s = node_offset(id);
            size_t l = xgidx->get_length(xgidx->get_handle(id));
            cairo_move_to(cr, s, y_pos);
            cairo_set_line_width(cr, 1);
            cairo_line_to(cr, s+l, y_pos);
            cairo_stroke(cr);
            // DNA
            if (show_dna) {
                cairo_text_extents_t te;
                cairo_select_font_face(cr, "Arial",
                                       CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
                cairo_set_font_size(cr, 1);
                string seq = xgidx->get_sequence(xgidx->get_handle(id));
                for (int i = 0; i < seq.size(); ++i) {
                    string c = string(1, seq[i]);
                    cairo_text_extents(cr, c.c_str(), &te);
                    cairo_move_to(cr, s+i+0.05, y_pos+2.1);
                    cairo_show_text(cr, c.c_str());
                }
            }
            // edges from
            xgidx->follow_edges(h, false, [&](const handle_t& o) {
                    id_t id2 = xgidx->get_id(o);
                    double s2 = node_offset(id2);
                    double x = s+l;
                    double y = y_pos;
                    int delta = s2 - x;
                    double w = pow(log(abs(delta)+1), 1.5);
                    int xdiff = (delta < 0 ? -w : w)/2;
                    int ydiff = w*2;
                    double x1 = x+xdiff, y1=y-ydiff,
                        x2 = s2-xdiff, y2=y-ydiff,
                        x3 = s2, y3 = y;
                    cairo_move_to(cr, x, y);
                    cairo_curve_to(cr, x1, y1, x2, y2, x3, y3);
                    cairo_set_line_width(cr, 0.25);
                    cairo_stroke(cr);
                    return true;
                });
        });
    y_pos += 4;
    if (show_paths) {
        xgidx->for_each_path_handle([&](const path_handle_t& path) {
                string path_name = xgidx->get_path_name(path);
                // write the path name
                {
                    set_hash_color(path_name);
                    cairo_text_extents_t te;
                    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
                    cairo_set_font_size(cr, 1);
                    cairo_text_extents(cr, path_name.c_str(), &te);
                    cairo_move_to(cr, left_border-(te.width+1), y_pos+0.35);
                    cairo_show_text(cr, path_name.c_str());
                }
                // determine counts
                map<id_t, int> counts;
                int max_count = 0;
                xgidx->for_each_step_in_path(path, [&](const step_handle_t& step) {
                        handle_t handle = xgidx->get_handle_of_step(step);
                        auto& c = counts[xgidx->get_id(handle)];
                        ++c;
                        max_count = max(c, max_count);
                    });
                for (auto& c : counts) {
                    //set_hash_color(path_name);
                    //for (auto& m : p.mapping()) {
                    // get the start and end position of the node relative to the image

                    id_t id = c.first;
                    double s = node_offset(id);
                    size_t l = xgidx->get_length(xgidx->get_handle(id));
                    for (int j = 0; j < (!show_cnv ? 1 : c.second); ++j) {
                        cairo_move_to(cr, s, y_pos+(j*2));
                        cairo_set_line_width(cr, 1);
                        cairo_line_to(cr, s+l, y_pos+(j*2));
                        cairo_stroke(cr);
                    }
                    if (!show_cnv && c.second > 1) {
                        cairo_text_extents_t te;
                        auto q = hash_to_rgb(path_name, 0.5);
                        cairo_set_source_rgb(cr, 1.0-get<0>(q), 1.0-get<1>(q), 1.0-get<2>(q));
                        cairo_select_font_face(cr, "Arial",
                                               CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
                        cairo_set_font_size(cr, 1);
                        stringstream ss;
                        ss << c.second << "X";
                        string label = ss.str();
                        cairo_text_extents(cr, label.c_str(), &te);
                        cairo_move_to(cr, s, y_pos+0.25);
                        cairo_show_text(cr, label.c_str());
                        set_hash_color(path_name); // undo color change
                    }
                }
                max_count = (!show_cnv ? 1 : max_count);
                y_pos += 2*max_count;
            });
    }
    y_pos += 2;
    // coverage maps
    for (int i = 0; i < packs->size(); ++i) {
        auto& pack = packs->at(i);
        auto& name = pack_names[i];
        set_hash_color(name);
        {
            cairo_text_extents_t te;
            cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
            cairo_set_font_size(cr, 1);
            cairo_text_extents(cr, name.c_str(), &te);
            cairo_move_to(cr, left_border-(te.width+1), y_pos);
            cairo_show_text(cr, name.c_str());
        }
        // for each position draw a bar representing coverage
        size_t max_coverage = 0;
        for (size_t j = 0; j < pack.coverage_size(); ++j) {
            max_coverage = max(max_coverage, pack.coverage_at_position(j));
        }
        // use max_coverage to normalize into 0/1
        for (size_t j = 0; j < pack.coverage_size(); ++j) {
            double c = (double)pack.coverage_at_position(j)/(double)max_coverage;
            double x_pos = j + left_border + nodes_before_offset(j) + 0.5;
            //cerr << "cov " << name << " " << j << " " << c << " " << x_pos << " " << y_pos << endl;
            // draw a line that tall
            cairo_move_to(cr, x_pos, y_pos);
            //cairo_set_line_width(cr, c*2);
            cairo_set_line_width(cr, 1);
            cairo_line_to(cr, x_pos, y_pos-c);
            cairo_stroke(cr);
        }
        y_pos += 2;
    }
}

void Viz::draw(void) {
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);
    draw_graph();
}

void Viz::close(void) {
    if (cr != nullptr && surface != nullptr) {
        if (output_png) {
            cairo_status_t status = cairo_surface_write_to_png(surface, outfile.c_str());
            switch (status) {
            case CAIRO_STATUS_SUCCESS:
                // No problem
                break;
            case CAIRO_STATUS_NO_MEMORY:
                throw std::runtime_error("There is not enough memory to save a " + std::to_string(image_width) + "x" + std::to_string(image_height) + " PNG!");
                break;
            case CAIRO_STATUS_SURFACE_TYPE_MISMATCH:
                throw std::runtime_error("The surface does not have pixel contents! We cannot save it as a PNG!");
                break;
            case CAIRO_STATUS_WRITE_ERROR:
                throw std::runtime_error("Encountered an I/O error when saving PNG");
                break;
            case CAIRO_STATUS_PNG_ERROR:
                throw std::runtime_error("libpng returned an error when saving PNG");
                break;
            default:
                throw std::runtime_error("Encountered an unrecognized Cairo error (" + std::string(cairo_status_to_string(status)) + ") when saving " + std::to_string(image_width) + "x" + std::to_string(image_height) + " PNG");
                break;
            }
        }
        cairo_destroy(cr);
        cairo_surface_destroy(surface);
        cr = nullptr;
        surface = nullptr;
    } else {
        // We should already be closed out
        assert(cr == nullptr);
        assert(surface == nullptr);
    }
}

}
