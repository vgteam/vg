#include <iostream>
#include <fstream>
#include "pb2json.h"
#include "vg.pb.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

int main(int argc,char *argv[])
{
    Node n;
    n.set_sequence("AAAACACGAAGGTGCTCACTGAAACATGGGAACCAAAGTTTCCCACAACATAAGGAGCAGAGTGAAACTGCAGAGGTACAATCATGGAATTCCAGAAAATG");
    n.set_name("14:20000000-20000100");
    n.set_id(1);
    ofstream of("test.vg");
    n.SerializeToOstream(&of);
	char *json2 = pb2json(n);
	cout<<json2<<endl;
	free(json2);
    return 0;
}
