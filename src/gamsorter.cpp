#include "gamsorter.hpp"
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

class AlnSorter{
    bool invert = false;
public:
  bool mycomparison(const bool& reverse=false)
  {invert = reverse;}
  bool operator() (const Alignment& a_one, const Alignment& a_two) const
  {

        // Unmapped case
        int xf_size = a_one.path().mapping_size();
        int xb_size = a_two.path().mapping_size();
        
        if (xf_size == 0 || xb_size == 0){
            bool both_zero = (xf_size == xb_size);
            if (both_zero){
                return false;
            }
            else{
                return (xf_size > xb_size) ? false : true;
            }
        }

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
};


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
        int xf_size = a_one.path().mapping_size();
        int xb_size = a_two.path().mapping_size();
        
        if (xf_size == 0 || xb_size == 0){
            bool both_zero = (xf_size == xb_size);
            if (both_zero){
                return false;
            }
            else{
                return (xf_size > xb_size) ? false : true;
            }
        }
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

void GAMSorter::write_temp(vector<Alignment>& alns)
{
    string t_name = "tmp_sort_" + std::to_string(tmp_filenames.size());
    tmp_filenames.push_back(t_name);
    ofstream t_file;
    t_file.open(t_name);
    stream::write_buffered(t_file, alns, 0);
    t_file.close();
}

void GAMSorter::dumb_sort(string gamfile)
{
    std::vector<Alignment> buf;
    buf.reserve(1000000);

    std::function<void(Alignment&)> presort = [&](Alignment &aln) {
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

void GAMSorter::stream_sort(string gamfile){

    std::vector<Alignment> buffer;
    buffer.reserve(1000);
    std::function<void(Alignment&)> firstpass_bufsort = [&](Alignment& a){
        #pragma omp critical
        buffer.push_back(a);
        if (buffer.size() == max_buf_size){
                sort(buffer);
                write_temp(buffer);
                buffer.clear();
        }
    };

    


    ifstream gammy;
    gammy.open(gamfile);
    stream::for_each(gammy, firstpass_bufsort);

    sort(buffer);
    write_temp(buffer);
    buffer.clear();



   
    int num_valid = tmp_filenames.size();
    vector<bool> still_valid(tmp_filenames.size());
    vector<uint64_t> curr_file_count (tmp_filenames.size());
    typedef priority_queue<Alignment,std::vector<Alignment>, AlnSorter> AlnHeap;
    

    for (int i = 0; i < tmp_filenames.size(); ++i){
        ifstream tfile;
        tfile.open(tmp_filenames[i]);
        tmp_files.emplace_back(&tfile);
    }


    

    vector<::google::protobuf::io::CodedInputStream*> protostreams;
    vector<::google::protobuf::io::GzipInputStream*> zipstreams;

    for (auto x : tmp_files){
        ::google::protobuf::io::ZeroCopyInputStream *raw_in =
          new ::google::protobuf::io::IstreamInputStream(x);
        ::google::protobuf::io::GzipInputStream *gzip_in =
          new ::google::protobuf::io::GzipInputStream(raw_in);
        ::google::protobuf::io::CodedInputStream *pstream =
          new ::google::protobuf::io::CodedInputStream(gzip_in);

        protostreams.push_back(pstream);
        zipstreams.push_back(gzip_in);
    }

    std::function<bool(int)> handle_count_line = [&](int tmp_index){
        uint64_t count;
        protostreams[tmp_index]->ReadVarint64((::google::protobuf::uint64*) &count);
        curr_file_count[tmp_index] = count;
        cerr << std::to_string(count) << endl;
        if (!count){
            return false;
        }
        return true;
    };

    for (int i = 0; i < protostreams.size(); i++){
        handle_count_line(i);
    }

    
    
    std::function<bool(int, Alignment&)> read_dat_buf = [&](int tmp_index, Alignment& ret){

        // Read our message of size count.
        std::string s;
        //for (uint64_t i = 0; i < count; ++i) {
            (curr_file_count[tmp_index])--;
            if (curr_file_count[tmp_index] <= 0){
                return handle_count_line(tmp_index);
            }
            uint32_t msgSize = 0;
            delete (protostreams[tmp_index]);
            (protostreams[tmp_index]) = new ::google::protobuf::io::CodedInputStream( (zipstreams[tmp_index]));
            // the messages are prefixed by their size
            (protostreams[tmp_index])->ReadVarint32(&msgSize);
            if ((msgSize > 0) &&
                (protostreams[tmp_index]->ReadString(&s, msgSize))) {
                ret.ParseFromString(s);
                return true;
            }
        //}
     };


    // Heap of size (tmp_files.size())
    AlnHeap lil_heap;
    for(int i = 0; i < tmp_files.size(); i++){
        Alignment a;
        read_dat_buf(i, a);
        lil_heap.push( a );
    }

    
    // Rotating index within the tmp_files array,
    // for accessing the ifstream and associated Proto readers/writers
    int tmp_file_index = 0;

    // This is just an output buffer.
    vector<Alignment> sorted_out_buf;
    sorted_out_buf.reserve(1000);

    // Our output stream
    string outname = gamfile + ".sorted.gam";
    ofstream ofile;
    ofile.open(outname);

    cerr << "in the loop" << endl;

    while(!lil_heap.empty()){

        Alignment a;
        bool valid = read_dat_buf(tmp_file_index, a);
        if (valid){
            lil_heap.push(a);
        }
        else{
            //handle_count_line(tmp_file_index);
        }
        

        
        sorted_out_buf.push_back(lil_heap.top());
        cerr << "writing " << lil_heap.top().name() << endl;
        if (sorted_out_buf.size() % 100 == 0){
            cerr << "Done " << sorted_out_buf.size() << " records" << endl;
        }
        lil_heap.pop();
        stream::write_buffered(ofile, sorted_out_buf, 1000);
        tmp_file_index = (tmp_file_index + 1) % tmp_files.size();

    }
    cerr << "out the loop" << endl;

    // write any remaining records in our buffer
    stream::write_buffered(ofile, sorted_out_buf, 0);

}

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
