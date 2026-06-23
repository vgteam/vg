#include "viz.hpp"
#include <regex>
#include <utility>

#include <dlfcn.h>

// To make Cairo an optional dependency, we include its headers but we wrap all
// its functions with wrappers that grab the functions from the runtime dynamic
// linker if possible and error otherwise. We need to do backflips to let the
// code be written as if this is not happening when we actually go to *use*
// Cairo.

#ifdef __APPLE__
#define LIBRARY_EXT "dylib"
#else
#define LIBRARY_EXT "so"
#endif

namespace vg {


#ifdef CAIRO_IS_OPTIONAL

std::atomic<void*> Viz::libcairo_handle = nullptr;

/**
 * Try to load Cairo if available and return the dlopen() handle to it. If not
 * possible, return nullptr.
 */
static void* load_cairo() {
    return dlopen("libcairo." LIBRARY_EXT, RTLD_NOW | RTLD_GLOBAL);
}

/**
 * Must be called before any Cairo functions are called.
 */
static void ensure_cairo_loaded() {
    if (Viz::libcairo_handle == nullptr) {
        void* cairo_handle = load_cairo();
        void* expected = nullptr;
        if (!Viz::libcairo_handle.compare_exchange_strong(expected, cairo_handle)) {
            dlclose(cairo_handle);
        }
    }
}

// We wrap each function symbol in Cairo that we use with an instantiation of
// this function for it.
#define MAKE_WRAPPER(fn_name) \
template<typename... Args> \
typename std::invoke_result<decltype(::fn_name), Args...>::type WRAPPED_ ## fn_name(Args... args) { \
    if (Viz::libcairo_handle == nullptr) { throw std::runtime_error("Cairo not available"); } \
    decltype(&::fn_name) resolved = (decltype(&::fn_name)) dlsym(Viz::libcairo_handle, #fn_name); \
    if (resolved == nullptr) { throw std::runtime_error(#fn_name " not available"); } \
    return (*resolved)(std::forward<Args>(args)...); \
}

// If any Cairo symbols are undefined at link time, add wrapping for them here.
// TODO: Figure out a way to avoid needing to redefine the function name for
// later code; defining a function of the same name in the vg namespace isn't
// it.
MAKE_WRAPPER(cairo_create)
#define cairo_create WRAPPED_cairo_create
MAKE_WRAPPER(cairo_curve_to)
#define cairo_curve_to WRAPPED_cairo_curve_to
MAKE_WRAPPER(cairo_destroy)
#define cairo_destroy WRAPPED_cairo_destroy
MAKE_WRAPPER(cairo_image_surface_create)
#define cairo_image_surface_create WRAPPED_cairo_image_surface_create
MAKE_WRAPPER(cairo_line_to)
#define cairo_line_to WRAPPED_cairo_line_to
MAKE_WRAPPER(cairo_move_to)
#define cairo_move_to WRAPPED_cairo_move_to
MAKE_WRAPPER(cairo_paint)
#define cairo_paint WRAPPED_cairo_paint
MAKE_WRAPPER(cairo_select_font_face)
#define cairo_select_font_face WRAPPED_cairo_select_font_face
MAKE_WRAPPER(cairo_set_font_size)
#define cairo_set_font_size WRAPPED_cairo_set_font_size
MAKE_WRAPPER(cairo_set_line_width)
#define cairo_set_line_width WRAPPED_cairo_set_line_width
MAKE_WRAPPER(cairo_set_source_rgb)
#define cairo_set_source_rgb WRAPPED_cairo_set_source_rgb
MAKE_WRAPPER(cairo_show_text)
#define cairo_show_text WRAPPED_cairo_show_text
MAKE_WRAPPER(cairo_status_to_string)
#define cairo_status_to_string WRAPPED_cairo_status_to_string
MAKE_WRAPPER(cairo_stroke)
#define cairo_stroke WRAPPED_cairo_stroke
MAKE_WRAPPER(cairo_surface_destroy)
#define cairo_surface_destroy WRAPPED_cairo_surface_destroy
MAKE_WRAPPER(cairo_surface_status)
#define cairo_surface_status WRAPPED_cairo_surface_status
MAKE_WRAPPER(cairo_surface_write_to_png)
#define cairo_surface_write_to_png WRAPPED_cairo_surface_write_to_png
MAKE_WRAPPER(cairo_svg_surface_create)
#define cairo_svg_surface_create WRAPPED_cairo_svg_surface_create
MAKE_WRAPPER(cairo_text_extents)
#define cairo_text_extents WRAPPED_cairo_text_extents

// Now Cairo code should work as normal.

#endif

Viz::Viz(PathHandleGraph* x, vector<Packer>* p, const vector<string>& n, const string& o, int w, int h, bool c, bool d, bool t) {
#ifdef CAIRO_IS_OPTIONAL
    ensure_cairo_loaded();    
#endif
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
        check_status(cairo_surface_status(surface), "creating " + std::to_string(image_width) + "x" + std::to_string(image_height) + " vector surface");
    } else {
        // Cairo can only handle a maximum size of 32,767 px on a side for a raster image.
        // See https://stackoverflow.com/a/24600155/402891
        if (image_width > 32767 || image_height > 32767) {
            // Complain to the user and stop.
            // TODO: Work out how to draw bigger images than Cairo will allow.
            std::cerr << "[Viz::init] error: Image would be " << image_width << "x" << image_height << " which is larger than the maximum 32,767 px on a side canvas that Cairo supports for raster images. Use an SVG output, try `odgi viz`, or reduce your input size!" << std::endl;
            exit(1);
        }
        surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, image_width, image_height);
        check_status(cairo_surface_status(surface), "creating " + std::to_string(image_width) + "x" + std::to_string(image_height) + " image surface");
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
    check_status(cairo_surface_status(surface), "creating zero-sized scratch surface");
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
            check_status(status, "saving " + std::to_string(image_width) + "x" + std::to_string(image_height) + " PNG");
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

void Viz::check_status(const cairo_status_t& status, const std::string& task) {
    switch (status) {
    case CAIRO_STATUS_SUCCESS:
        // No problem
        break;
    case CAIRO_STATUS_NO_MEMORY:
        throw std::runtime_error("Ran out of memory while " + task);
        break;
    case CAIRO_STATUS_INVALID_SIZE:
        throw std::runtime_error("The input was too big when " + task);
        break;
    case CAIRO_STATUS_SURFACE_TYPE_MISMATCH:
        throw std::runtime_error("The surface does not have pixel contents! We cannot attempt " + task);
        break;
    case CAIRO_STATUS_WRITE_ERROR:
        throw std::runtime_error("Encountered an I/O error while " + task);
        break;
    case CAIRO_STATUS_PNG_ERROR:
        throw std::runtime_error("libpng returned an error while " + task);
        break;
    default:
        throw std::runtime_error("Encountered an unrecognized Cairo error (" + std::string(cairo_status_to_string(status)) + ") while " + task);
        break;
    }
}

}
