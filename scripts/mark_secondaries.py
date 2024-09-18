#!/usr/bin/python3
# mark_secondaries.py: Mark all but the first alignment with a given name as secondary
"""
Mark duplicate alignments for a given read name as secondary. Useful for GraphAligner output which does not mark its secondaries. Assumes that the first alignment is the primary alignment, ignoring score.

  vg view -a unmarked.gam mark_secondaries.py | vg view -JaG - > marked.gam

"""
import sys
import json


def filter_json_gam(infile):
    """
    process gam json made with vg view -a my.gam
    """

    seen_names = set()
    
    for line in infile:
        gam = json.loads(line)
        
        if gam['name'] in seen_names:
            gam['is_secondary'] = True
        else:
            gam['is_secondary'] = False
            seen_names.add(gam['name'])

        print(json.dumps(gam))

def main():
    """
    Main entry point for the program.
    """
    filter_json_gam(sys.stdin)

if __name__ == "__main__" :
    main()
        
