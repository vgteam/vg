import os
import sys
import shutil


def make_source_me():
    var_dict = {
            "LIBRARY_PATH":"lib",
            "LD_LIBRARY_PATH":"lib",
            "LD_INCLUDE_PATH":"include",
            "C_INCLUDE_PATH":"include",
            "CPLUS_INCLUDE_PATH":"include",
            "INCLUDE_PATH":"include",
            "PATH":"bin"                                                                   }

    my_wd = os.getcwd()
    out_file = "source_me.sh"

    with open(out_file, "w") as ofi:
        for i in var_dict:
            o_line  = "export" + " " + i + "=" + my_wd + "/" + var_dict[i] + ":$" + i + "\n"
            ofi.write(o_line)

def set_compiler():
    compiler_dict = {
        "CC": "gcc",
        "CXX": "g++"
        }

    out_file = "source_me.sh"
    with open(out_file, "a") as ofi:
        for i in compiler_dict:
            o_line = "export" + " " + i + "=" + "$(which " + compiler_dict[i] + ")\n"
            ofi.write(o_line)

def check_deps():
    return

## pulled these out with this line: for i in `cat source_me.sh | cut -f 2 -d " " | cut -f 1 -d "="`; do echo "\"\$$i\":\"\","; done
if __name__ == "__main__":
    make_source_me() 
    set_compiler()
## "$CC":"",
            #"$CXX":"",  
