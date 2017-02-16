#include "GAMSorter.hpp"
/*
*  GAMSorter: sort a gam by position and offset
*  dumbly store unmapped reads at the end.
*/

using namespace std;
using namespace vg;

/**
vector<Alignment> GAMSorter::merge(vector<vector<Alignment> > a){
    vector<Alignment> ret;
    vector<int> a_tracker(a.size());
    


    return ret;
}

void GAMSorter::merge(map<int, vector<Alignment> > m, map<int, int> split_to_sz){
    
}

void GAMSorter::merge(vector<string> sorted_tmp_filenames){

}

**/

struct custom_pos_sort_key
{
    bool operator()(const Position lhs, const Position rhs)
    {
        if (lhs.node_id() == rhs.node_id())
        {
            return lhs.offset() < rhs.offset();
        }
        else
        {
            return lhs.node_id() < rhs.node_id();
        }
    }
} possortkey;

struct custom_aln_sort_key
{
    bool operator()(const Alignment a_one, const Alignment a_two)
    {
        Position x_front = a_one.path().mapping(0).position();
        Position x_back = a_one.path().mapping(a_one.path().mapping_size() - 1).position();

        Position lhs;
        Position rhs;
        if (x_front.node_id() < x_back.node_id())
        {
            lhs = x_front;
        }
        else
        {
            lhs = x_back;
        }
        Position y_front = a_two.path().mapping(0).position();
        Position y_back = a_two.path().mapping(a_two.path().mapping_size() - 1).position();
        if (y_front.node_id() < y_back.node_id())
        {
            rhs = y_front;
        }
        else
        {
            rhs = y_back;
        }

        if (lhs.node_id() == rhs.node_id())
        {
            return lhs.offset() < rhs.offset();
        }
        else
        {
            return lhs.node_id() < rhs.node_id();
        }
    }
} alnsortkey;

void GAMSorter::sort(vector<Alignment> &alns)
{
    std::sort(alns.begin(), alns.end(), alnsortkey);
}

void GAMSorter::paired_sort(string gamfile)
{

    std::unordered_map<string, Alignment> pair_seconds;
    vector<Alignment> firsts;
    firsts.reserve(100000);
    std::function<void(Alignment &, Alignment &)> sort_the_gam = [&](Alignment &alpha, Alignment &beta) {
        bool first_min = min_aln_first(alpha, beta);
        if (!alpha.fragment_next().name().empty())
        {
            if (first_min)
            {
                firsts.push_back(alpha);
                pair_seconds[beta.name()] = beta;
            }
            else
            {
                firsts.push_back(beta);
                pair_seconds[alpha.name()] = alpha;
            }
        }
        else
        {
            firsts.push_back(alpha);
        }

    };

    ifstream gammy;
    gammy.open(gamfile);
    stream::for_each_interleaved_pair_parallel(gammy, sort_the_gam);

    sort(firsts);
    vector<Alignment> ret;
    for (auto x : firsts)
    {
        // write first read.
        ret.push_back(x);
        // grab second read and write it.
        if (!x.fragment_next().name().empty())
            ret.push_back(pair_seconds[x.fragment_next().name()]);
    }
}

void GAMSorter::write_temp(vector<Alignment> &alns)
{
    string t_name = "tmp_sort_" + tmp_filenames.size();
    tmp_filenames.push_back(t_name);
    ofstream t_file;
    t_file.open(t_name);
    stream::write_buffered(t_file, alns, 1000);
}

void GAMSorter::dumb_sort(string gamfile)
{
    std::vector<Alignment> buf;
    buf.reserve(1000000);

    std::function<void(Alignment &)> presort = [&](Alignment &aln) {

        buf.push_back(aln);
    };
    ifstream gammy;
    gammy.open(gamfile);
    stream::for_each(gammy, presort);

    std::sort(buf.begin(), buf.end(), alnsortkey);

    std::function<void(uint64_t)> x_buf_func = [&](uint64_t) {

    };

    ofstream outfi;
    outfi.open(gamfile + ".sorted.gam");
    stream::write_buffered(outfi, buf, buf.size());
}

// void GAMSorter::stream_sort(string gamfile){

//     std::vector<Alignment&> buffer;
//     buffer.reserve(1000);
//     std::function<void(Alignment&)> firstpass_bufsort = [&]{
//         #pragma omp critical
//         buffer.push_back(a);
//         if (buffer.size() == max_buf_size){
//                 sort(buffer);
//                 write_temp(buffer);
//                 buffer.clear();
//         }

//     };

//     ifstream gammy;
//     gammy.open(gamfile);
//     stream::for_each_interleaved_pair_parallel(gammy, firstpass_bufsort);

//     vector<Alignment> file_tops(tmp_files.size());
//     std::function<int(vector<Alignment>)> min_index = [&](vector<Alignment> alns){
//         int min_ind = 0;
//         for (int i = 0; i < alns.size(); ++i){
//             if (alns[i] < alns[min_ind]){
//                 min_ind = i;
//             }
//         }
//         return i;
//     };

//     std::function<void(int, vector<Alignment>&)> handle_lowest_tmpfile = [&](int min_ind,
//                                                                             ofstream ofi,
//                                                                             vector<Alignment> buf,
//                                                                             vector<Alignment>& alns){
//         //ofi.writeline(alns[i]);
//         //alns[i] = tmp_files[min_ind].getline();
//         if (tmp_files[min_ind].eof())){
//             tmp_files.erase(tmp_files.begin() + min_ind);
//             alns.erase(alns.begin() + min_ind);
//         }
//         if (buf.size() >= max_buf_size){
//             stream::write_buffered(ofi, buf, 1000);
//         }
//         };

//     /**
//     * Open all our files and maintain a buffer of size N(tmp_files)
//     * at each iteration, perform a merge and write to disk
//     */
//     for (auto tmpfi : tmp_file_names){
//             // Open all our tmp_files and load 1 record into
//             // our buffer
//             tmp_files.emplace_back(ifstream(tmpfi));

//     }

//     ofstream outfi;
//     outfi.open("x.txt");
//     int min_ind = 0;
//     vector<Alignment> obuf;
//     while(file_tops.size() > 0){
//         min_ind = min_index(file_tops);
//         handle_lowest_tmpfile(min_ind, outfi, obuf, file_tops);
//     }

// }

void GAMSorter::write_index(string gamfile, string outfile, bool isSorted)
{

    map<int64_t, vector<string>> node_to_alns;

    std::function<void(Alignment &)> node_grabber = [&](Alignment &aln) {
        vector<int64_t> ret;
        Path p = aln.path();
        for (int i = 0; i < p.mapping_size(); i++)
        {
            Mapping m = p.mapping(i);
            ret.push_back(m.position().node_id());
        }
        for (auto x : ret)
        {
            node_to_alns[x].push_back(aln.name());
        }
    };
}

bool GAMSorter::min_aln_first(Alignment &a, Alignment &b)
{

    if (less_than(get_min_position(a), get_min_position(b)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

Position GAMSorter::get_min_position(Alignment a)
{
    return get_min_position(a.path());
}

Position GAMSorter::get_min_position(Path pat)
{
    Position p = pat.mapping(0).position();
    Position p_prime = pat.mapping(pat.mapping_size() - 1).position();
    return less_than(p, p_prime) ? p : p_prime;
}

bool GAMSorter::equal_to(Position a, Position b)
{
    if (a.node_id() == b.node_id() &&
        a.offset() == b.offset())
    {
        return true;
    }
    return false;
}

bool GAMSorter::less_than(Position a, Position b)
{
    if (a.node_id() > b.node_id())
    {
        return false;
    }
    else
    {
        return (a.node_id() == b.node_id()) ? (a.offset() < b.offset()) : true;
    }
}

bool GAMSorter::greater_than(Position a, Position b)
{
    if (a.node_id() < b.node_id())
    {
        return false;
    }
    else
    {
        return (a.node_id() == b.node_id()) ? (a.offset() > b.offset()) : true;
    }
}