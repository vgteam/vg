#ifndef GAMSORT_H
#define GAMSORT_H

#include "vg.pb.h"
#include "stream.hpp"
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <unordered_map>

/**
 Gam sort
*/
using namespace std;
namespace vg
{



class GAMSorter
{

  
  public:
    // vector<Alignment> merge(vector<vector<Alignment>> a);

    // void merge(map<int, vector<Alignment>> m, map<int, int> split_to_sz);

    void sort(vector<Alignment>& alns);

    void paired_sort(string gamfile);

    void write_temp(vector<Alignment>& alns);

    // void stream_sort(string gamfile);

    void dumb_sort(string gamfile);

    // vector<Alignment> split(vector<Alignment> a, int s);

    bool min_aln_first(Alignment& a, Alignment& b);

    Position get_min_position(Alignment a);

    Position get_min_position(Path p);

    void write_index(string gamfile, string outfile, bool isSorted = false);

    bool equal_to(Position a, Position b);

    bool less_than(Position a, Position b);

    bool greater_than(Position a, Position b);

  private:
    map<int, int> split_to_split_size;
    vector<string> tmp_filenames;
    vector<ifstream> tmp_files;
    int max_buf_size = 20000;
    /**
    * We want to keep pairs together, with the lowest-coordinate pair coming first.
    * If one read is unmapped, it follows its partner in the sorted GAM file.
    * Pairs with two unmapped reads are considered the highest-mapped in the graph, i.e.
    * they come at the end of a sorted GAM file.

    * Since we can't hash alignments, we store both mates in the unordered map by name.
    */
    unordered_map<string, Alignment> pairs;
    unordered_map<string, pair<Alignment, Alignment> > paired_pairs;
   
};
}
#endif
