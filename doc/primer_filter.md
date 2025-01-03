
`vg primers` can be used to filter PRC primers based on properties of the pangenome such as whether there are variations in the primers and the possible lengths of the PRC product.
`vg primers` takes pangenome indexes and the output of `primer3` as input and outputs a `.tsv` file of the input primers and properties from the pangenome.

# Get primers with `primer3`

The input to `vg primers` is the output of the command line version of [`primer3`](https://github.com/primer3-org/primer3).

`primer3` requires a config file formatted like:
```
SEQUENCE_ID=CHM13#0#chr17|BRCA1P1|exon_1|44026826
SEQUENCE_TEMPLATE=CATGT...
PRIMER_NUM_RETURN=10
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
=
```

The `SEQUENCE_ID` field must be formatted correctly for `vg primers` to find the correct location of the primers in the pangenome.
There are four fields in `SEQUENCE_ID` separated by `|`. 
They are the reference path name in the graph, the name of the gene, the exon or intron, and the offset of the sequence in the path.
The names of the reference paths in the graph can be found with `vg paths -L -R`.
The `SEQUENCE_TEMPLATE` is the nucleotide sequence the primers are found from.

# Filtering primers with `vg primers`

`vg primers` requires the following indexes of the pangenome:

- the *xg* index created with `vg index -x`
- the *distance*/*snarl* index created with `vg index -j`
- the [*r-index*](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand) created with `vg gbwt -r`
- the [*gbz*](https://github.com/vgteam/vg/wiki/VG-GBWT-Subcommand) created with `vg gbwt`

# Interpreting output of `vg primers`

`vg primers` outputs a tsv file with the following fields for each primer:

| field      | definition | description |
| ---------- | ---------- | ----------- |
| chrom      | chromosome | reference path name, the first field in `SEQUENCE_ID` |
| tplfeat    | template feature | the second and third fields in `SEQUENCE_ID` |
| tplpos     | template position | offset along the reference path, the fourth field in `SEQUENCE_ID` |
| lpseq      | left primer sequence | the nucleotide sequence of the left primer |
| rpseq      | right primer sequence | the nucleotide sequence of the right primer |
| lppostpl   | left primer position template | the offset of the left primer in the template sequence |
| rppostmp   | right primer position template | the offset of the right primer in the template sequence |
| lpposchrom | left primer position chromosome | the offset of the left primer in the reference path |
| rpposchrom | right primer position chromosome | the offset of the right primer in the reference path |
| pnid       | left primer mapped node ids | the node ids that the left primer overlaps in the graph |
| rpnid      | right primer mapped node ids | the node ids that the right primer overlaps in the graph |
| lplen      | left primer length | the length in nucleotides of the left primer |
| rplen      | right primer length | the length in nucleotides of the right primer |
| linsize    | linear product size | the length of the product in the linear genome (including primer lengths) | 
| minsize      | minimum product size | the minimum length of the product according to the pangenome |
| maxsize      | maximum product size | the maximum length of the product according to the pangenome |
| varlevel     | variation level | a measure of variation in the primers (the number of unique haplotypes / the total number of haplotypes) |

