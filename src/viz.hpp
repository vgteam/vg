#ifndef VG_VIZ_HPP_INCLUDED
#define VG_VIZ_HPP_INCLUDED

#include <vector>
#include <cairo.h>
#include <cairo-svg.h>
#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>
#include "hash_map.hpp"
#include "packer.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "vg/io/json2pb.h"

namespace vg {

using namespace std;

class Viz {
public:
    Viz(void) { }
    ~Viz(void) { close(); }
    Viz(PathHandleGraph* x, vector<Packer>* p, const vector<string>& n, const string& o, int w, int h, bool c, bool d, bool t);
    void init(PathHandleGraph* x, vector<Packer>* p, const vector<string>& n, const string& o, int w, int h, bool c, bool d, bool t);
    void draw(void);
    void draw_graph(void);
    void close(void);
private:
    double node_offset(id_t id);
    double nodes_before_offset(size_t pos);
    uint64_t id_to_rank(nid_t id);
    void set_hash_color(const string& str);
    void compute_borders_and_dimensions(void);
    void check_status(const cairo_status_t& status, const std::string& task = "calling Cairo");
    unordered_map<nid_t, uint64_t> id_rank_map;
    PathHandleGraph* xgidx = nullptr;
    vector<Packer>* packs = nullptr;
    vector<string> pack_names;
    string outfile;
    cairo_surface_t *surface = nullptr;
	cairo_t *cr = nullptr;
    bool output_png = false;
    bool output_svg = false;
    bool show_cnv = true;
    bool show_dna = true;
    bool show_paths = true;
    int image_width = 0;
    int image_height = 0;
    int left_border = 0;
    int top_border = 0;
    size_t seq_length = 0;
};

tuple<double, double, double> hash_to_rgb(const string& str, double min_sum);

}

#endif
