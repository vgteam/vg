# The Graph Alignment Format (GAF)

The latest version of this document is hosted at: <https://github.com/vgteam/vg/blob/master/doc/static/GAF.md>

This document describes **version 1.0** of the [vg](https://github.com/vgteam/vg) interpretation of the Graph Alignment Format (GAF).
It is a superset of a subset of the [original GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md).
That format in turn is a superset of the [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).
Sequence names and optional fields follow conventions set in the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).
Difference strings are defined in the [minimap2 man page](https://lh3.github.io/minimap2/minimap2.html).
Paths are represented as [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) walks.
Reference graphs may have [pggname](https://github.com/jltsiren/pggname) stable names.

## Overview

GAF is a tab-delimited file format for sequence alignments to bidirected sequence graphs.
The file is encoded in UTF-8.
Unless otherwise specified, all fields are restricted to 7-bit US-ASCII.

Each file consists of a number of header lines followed by a number of alignment lines.
Each line can be split into a number of fields separated by TAB (`\t`) characters.

## Typed fields

Typed fields are stored in the SAM-style `TAG:TYPE:VALUE` format.
The tag is a two-character string matching `[A-Za-z][A-Za-z0-9]`.

The following types are currently supported:

|Type|Description|
|:--:|:----------|
|`A` |Printable character in `[!-~]`|
|`Z` |String of printable characters and spaces (`[ !-~]*`)|
|`i` |Signed 64-bit integer|
|`f` |Double-precision floating point number|
|`b` |Boolean value, with `1` for true and `0` for false|

## Header lines

**Since vg 1.70.0**

Header lines are optional, and they must all appear before the first alignment line.
The first field of each header line is a three-character tag matching `@[A-Za-z][A-Za-z0-9]`.

**Example:**
```txt
@HD	VN:Z:1.0
@RN	7f4b28c71ceb808aebd8b8e9fe85e79d0d208ee263ffe9fcdef5ade20534ceb5
@SG	7f4b28c71ceb808aebd8b8e9fe85e79d0d208ee263ffe9fcdef5ade20534ceb5	e10f3b362d8a4273059d9aea38a78bd71913418c3f3c9a2b5ea44e86de2c1181
@TL	e10f3b362d8a4273059d9aea38a78bd71913418c3f3c9a2b5ea44e86de2c1181	1f133f116e8dd98fc07a647a8954038c2bcf07a45759ba94718471fe34ed7a7c
```

### File headers

File headers start with tag `@HD`.
They may contain any number of optional typed fields.
The following optional fields are known.

|Tag |Type|Description|
|:--:|:--:|:----------|
|`VN`|`Z` |Version number (e.g. `1.0`; only one allowed in the file)|

### Reference name

**Since vg 1.71.0**

The graph the sequences were aligned to can be identified using a reference name line.
A reference name line starts with tag `@RN` and contains the pggname (SHA-256 hash of the canonical GFA representation) of the graph as the second field.
There may be optional typed fields.
There can be only one reference name line in a file.

### Graph relationships

**Since vg 1.71.0**

If the sequences were aligned to graph A, which is a subgraph of B, graph B is also a valid reference for the alignments.
If there is a known coordinate translation from graph B to graph C, graph C can also be used as a reference after translating the coordinates.
Subgraph (`@SG`) and translation (`@TL`) lines can used to describe such relationships between reference graphs.

A subgraph line contains the pggname of the subgraph as the second field and the name of the supergraph as the third field.
A translation line contains the name of the source graph as the second field and the name of the destination graph as the third field.
Both line types may contain optional typed fields, and there may be any number of such lines.

## Alignment lines

Each alignment line has 12 mandatory fields.
Missing values in fields 3 to 11 are indicated by character `*`.

|Field|Type|Description|
|----:|:--:|:----------|
|1    |`Z` |Query sequence name|
|2    |`i` |Query sequence length|
|3    |`i` |Query start (0-based; closed)|
|4    |`i` |Query end (0-based; open)|
|5    |`A` |Strand relative to the path; always `+`|
|6    |`Z` |Target path represented as a GFA walk|
|7    |`i` |Target path length|
|8    |`i` |Start position on the target path (0-based; closed)|
|9    |`i` |End position on the target path (0-based; open)|
|10   |`i` |Number of matches|
|11   |`i` |Number of matches, mismatches, insertions, and deletions|
|12   |`i` |Mapping quality (0-255; 255 for missing)|

**Example:**
```txt
read1 	6 	0 	6 	+ 	>2>3>4 	12 	2 	8 	6 	6 	60 	cs:Z::6
read2 	7 	0 	7 	+ 	>2>5>6 	11 	1 	8 	7 	7 	60 	cs:Z::7
read3 	7 	0 	7 	* 	* 	* 	* 	* 	* 	* 	255	cs:Z:+GATTACA
```

### Query sequence name

Query sequence names must follow SAM conventions.
A name may contain any printable ASCII characters in the range `[!-~]`, except `@`.
This allows distinguishing header lines from alignment lines.

### Target path

This version of GAF does not allow specifying the target path using stable rGFA coordinates or nodes (GFA segments) with string names.
Nodes must have positive integer identifiers.
Node identifier `0` cannot be used, as many graph implementations reserve it for technical purposes.

### Optional fields

Optional fields are SAM-style typed fields.
No tag can appear more than once on the same line, and the order of the optional fields does not matter.

### Difference string

Difference strings represent an edit script that transforms the given interval of the target path to the given interval of the query sequence.
They are stored as an optional field `cs` of type `Z`.
We support a subset of the operations defined for minimap2 difference strings.

|Operation|Regex   |Description|
|:-------:|:------:|:----------|
|`:`      |`[0-9]+`|Number of matching bases|
|`*`      |`[ACGTN][ACGTN]`|Mismatch as (target base, query base)|
|`+`      |`[ACGTN]+`|Insertion as the unaligned query bases|
|`-`      |`[ACGTN]+`|Deletion as the unaligned target bases|

### Other defined optional fields

|Tag |Type|Description|
|:--:|:--:|:----------|
|`AS`|`i` |Alignment score|
|`bq`|`Z` |Base quality string; must have the same length as the query sequence|
|`fn`|`Z` |Name of the next fragment (for paired alignments; cannot be used with `fp`)|
|`fp`|`Z` |Name of the previous fragment (for paired alignments; cannot be used with `fn`)|
|`pd`|`b` |This alignment and its pair (specified by `fn` or `fp`) are properly paired|
|`fi`|`i` |Fragment identifier for a fragmented alignment (see below)|

## Conventions

### Header lines

**Since vg 1.70.0**

The first line of a GAF file is a file header (`@HD`) with the version number (`VZ:Z`) tag.
Any additional file header lines follow.
Other types of header lines are after file header lines.

### Primary alignments

A primary alignment represents an alignment of the entire query sequence to a non-empty interval of a target path.
Query start (field 3) must be `0` and query end (field 4) must have the same value as query sequence length (field 2).
A difference string must be present to allow recovering the entire query sequence.

### Unaligned sequences

**Since vg 1.70.0**

An unaligned sequence is represented as an alignment of the entire query sequence to a missing interval of a missing target path.
Query start (field 3) must be `0` and query end (field 4) must have the same value as query sequence length (field 2).
A difference string must be present, with the entire query sequence as a single insertion, to allow recovering the sequence.

### Fragmented alignments

A fragmented alignment is a single alignment represented as number of alignment lines (e.g. corresponding to subpaths that are within a specific subgraph).
The fragments (alignment lines) correspond to non-overlapping intervals of the underlying alignment.
Each fragment represents an alignment of a non-empty query interval to a non-empty interval of a non-empty target path.

Query interval (fields 3 and 4), target path (fields 6 to 9), and the difference string must be specific to each fragment.
Alignment statistics (fields 10 to 12) may be inherited from the underlying alignment or be specific to each fragment.
Fragments are identified by fragment indexes starting from `1`, stored as an optional field `fi` of type `i`.
