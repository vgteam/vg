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

   struct custom_pos_sort_key{
      bool operator() (const Position lhs, const Position rhs){
        if (lhs.node_id() == rhs.node_id()){
          return lhs.offset() < rhs.offset();
        }
        else{
          return lhs.node_id() < rhs.node_id();
        }
      }
    } possortkey;

    struct custom_aln_sort_key{
      bool operator() (const Alignment a_one, const Alignment a_two){
        Position x_front = a_one.path().mapping(0).position();
        Position x_back = a_one.path().mapping( a_one.path().mapping_size() - 1 ).position();

        Position lhs;
        Position rhs;
        if (x_front.node_id() < x_back.node_id()){
          lhs = x_front;
        }
        else{
          lhs = x_back;
        }
        Position y_front = a_two.path().mapping(0).position();
        Position y_back = a_two.path().mapping( a_two.path().mapping_size() - 1 ).position();
        if (y_front.node_id() < y_back.node_id()){
          rhs = y_front;
        }
        else{
          rhs = y_back;
        }




        if (lhs.node_id() == rhs.node_id()){
          return lhs.offset() < rhs.offset();
        }
        else{
          return lhs.node_id() < rhs.node_id();
        }
      }
    } alnsortkey;


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

   // pair<Alignment, Alignment> min_aln_first(Alignment a, Alignment b);

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
