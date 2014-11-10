#include <iostream>
#include <fstream>
#include "pb2json.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"
#include "vg.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

int main(int argc,char *argv[])
{

    vcf::VariantCallFile variantFile;

    string vcffilename = argv[1];
    variantFile.open(vcffilename);
    if (!variantFile.is_open()) {
        return 1;
    }

    string reffilename = argv[2];
    long int id = 0;

    vcf::Variant var(variantFile);
    while (variantFile.getNextVariant(var)) {
    }

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
