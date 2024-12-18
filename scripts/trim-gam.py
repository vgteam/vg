#!/usr/bin/env python3
"""
trim-gam.py: trim GAM records from JSON standard input to contain only mappings to nodes in the given range.

Range end is exclusive.

usage: vg view -aj input.gam | trim-gam.py 123 456 | vg view -JGa - >output.gam
"""

import sys
import json
import base64

def mapping_to_length(mapping: dict) -> int:
    """
    Get the total read bases involved in a mapping.
    """
    return sum(edit.get("to_length", 0) for edit in mapping.get("edit", []))

def trim_record(gam_record: dict, from_id: int, to_id: int) -> dict:
    """
    Trim down a single GAM record.
    """
    
    # Grab the sequence
    sequence = gam_record.get("sequence", "")
    # Decode the qualities. If empty, the result is empty
    quality = base64.b64decode(gam_record.get("quality", ""))

    # Collect together the sequence and quality pieces and the mappings we are keeping.
    sequence_parts: list[str] = []
    quality_parts: list[bytes] = []
    path_mappings: list[dict] = []

    # Start at read base 0
    read_cursor = 0

    for mapping in gam_record.get("path", {}).get("mapping", []):
        # Find the node each mapping maps to
        mapped_id = int(mapping.get("position", {}).get("node_id", "0"))
        # And the number of read bases used
        to_length = mapping_to_length(mapping)
        if mapped_id >= from_id and mapped_id < to_id:
            # We want to keep this mapping
            path_mappings.append(mapping)
            # And its read sequence
            sequence_parts.append(sequence[read_cursor:read_cursor + to_length])
            # And its quality. If there's no quality this is an empty string.
            quality_parts.append(quality[read_cursor:read_cursor + to_length])
        # Advance the cursor
        read_cursor += to_length
    
    # Copy the record
    new_record = dict(gam_record)

    # Install the trimmed pieces
    new_record["sequence"] = "".join(sequence_parts)
    new_record["quality"] = base64.b64encode(b"".join(quality_parts)).decode("utf-8")
    new_record.setdefault("path", {})["mapping"] = path_mappings

    return new_record

if __name__ == "__main__":

    if len(sys.argv) != 3:
        sys.stderr.write(__doc__)
        sys.stderr.write("\n")
        sys.exit(1)

    from_id = int(sys.argv[1])
    to_id = int(sys.argv[2])

    for line in sys.stdin:
        sys.stdout.write(json.dumps(trim_record(json.loads(line), from_id, to_id)))
        sys.stdout.write("\n")

