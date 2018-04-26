#include "viz.hpp"
#include <regex>

namespace vg {

Viz::Viz(xg::XG* x, vector<Packer>* p, const string& o, int w, int h, bool c) {
    init(x, p, o, w, h, c);
}

void Viz::init(xg::XG* x, vector<Packer>* p, const string& o, int w, int h, bool c) {
    xgidx = x;
    packs = p;
    outfile = o;
    left_border = 32;
    show_cnv = c;
    image_height = h;
    image_width = (w ? w : xgidx->seq_length + xgidx->node_count + left_border*2);
    top_border = image_height/2;
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

double Viz::node_offset(id_t id) {
    // provides 1 extra pixel of space per node plus border
    return xgidx->node_start(id) + xgidx->id_to_rank(id) + left_border;
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
    xgidx->for_each_handle([&](const handle_t& h) {
            // get the start and end position of the node relative to the image
            id_t id = xgidx->get_id(h);
            double s = node_offset(id);
            size_t l = xgidx->node_length(id);
            cairo_move_to(cr, s, top_border);
            cairo_set_line_width(cr, 1);
            cairo_line_to(cr, s+l, top_border);
            cairo_stroke(cr);
            /*
            // tick marks
            cairo_move_to(cr, s, 32);
            cairo_set_line_width(cr, 0.3);
            cairo_line_to(cr, s, 32+2);
            cairo_line_to(cr, s, 32-2);
            cairo_stroke(cr);
            */
            // edges from
            xgidx->follow_edges(h, false, [&](const handle_t& o) {
                    id_t id2 = xgidx->get_id(o);
                    double s2 = node_offset(id2);
                    double x = s+l;
                    double y = top_border;
                    int delta = s2 - x;
                    double w = pow(log(abs(delta)+1), 1.5);
                    int xdiff = (delta < 0 ? -w : w)/2;
                    int ydiff = w*2;
                    double x1 = x+xdiff, y1=y-ydiff,
                        x2 = s2-xdiff, y2=y-ydiff,
                        x3 = s2, y3 = y;
                    cairo_move_to(cr, x, y);
                    cairo_curve_to(cr, x1, y1, x2, y2, x3, y3);
                    cairo_set_line_width(cr, 0.5);
                    cairo_stroke(cr);
                    return true;
                });
        });
    int y_pos = top_border + 4;
    for (size_t i = 1; i <= xgidx->path_count; ++i) {
        string path_name = xgidx->path_name(i);
        set_hash_color(path_name);
        Path p = xgidx->path(path_name);
        // determine counts
        map<id_t, int> counts;
        int max_count = 0;
        for (auto& m : p.mapping()) {
            auto& c = counts[m.position().node_id()];
            ++c;
            max_count = max(c, max_count);
        }
        for (auto& c : counts) {
            //for (auto& m : p.mapping()) {
            // get the start and end position of the node relative to the image
            id_t id = c.first;//m.position().node_id();
            double s = node_offset(id);
            size_t l = xgidx->node_length(id);
            int node_copy_number = (!show_cnv ? 1 : c.second);
            for (int j = 0; j < c.second; ++j) {
                cairo_move_to(cr, s, y_pos+(j*2));
                cairo_set_line_width(cr, 1);
                cairo_line_to(cr, s+l, y_pos+(j*2));
                cairo_stroke(cr);
            }
        }
        max_count = (!show_cnv ? 1 : max_count);
        y_pos += 2*max_count;
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
            cairo_surface_write_to_png(surface, outfile.c_str());
        }
        cairo_destroy(cr);
        cairo_surface_destroy(surface);
        cr = nullptr;
        surface = nullptr;
    }
}

}
