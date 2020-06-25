import sys
import argparse


class FQRecord:
    def __init__(self):
        self.seq = ""
        self.qual = ""
        self.anno = ""
        self.name = ""

    def string(self):
        ret = "@" + self.name.strip() + "\n" + \
        self.seq + "\n" + \
        "+" + self.anno + "\n" + \
        self.qual
        return ret


"""Takes in a pseudo fasta and makes a
pseudo fastq from it"""
def make_fastq(record, name, fake_qual, seq="", anno=""):
    if record.name == "":
        record.name = name
    
    record.qual = "".join([fake_qual for i in range(0, len(record.seq))])

    return record

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="infile", type=str, required=True)
    parser.add_argument("-n", dest="namebase", type=str, required=True)
    parser.add_argument("-q", dest="qual", type=str, required=True)
    return parser.parse_args()
    

if __name__ == "__main__":
    ## Add the sample name to the front of the
    ## read line.
    args = parse_args()
    readfile = args.infile
    name_base = args.namebase
    fake_qual = args.qual
    count = 0
    with open(readfile, "r") as fi:
        for line in fi:
            record = FQRecord()
            record.seq = line.strip()
            record.name = name_base + "_" + str(count)
            record = make_fastq(record, "", fake_qual)
            print(record.string())
            count += 1
