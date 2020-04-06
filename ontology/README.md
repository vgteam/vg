# VG in RDF

## Conceptual model

`Node`s, `Path`s and `Step`s, are the three core parts of any VG graph in RDF.
A `Node` in the VG RDF corersponds directly to the Node concept in the VG protobuf serialization.
`Paths` are a number of `Step`s that represent a sequence of Node visits that generate a linear biological sequence.
Each `Step` connects a `Node` into a `Path`



## Annotations on a VG

Path are annotations on the graph, some paths represent a reference genome, others are the genome of a patient.
Yet some represent a gene, regulatory element or any other feature that we can annotate on any biological sequence.

For example lets say we have a small VG graph:

```turtle
node:100 a vg:Node ; vg:linksForwardToForward node:101 ; rdf:value "ACT" .
node:101 a vg:Node ; vg:linksForwardToForward node:102 : rdf:value "TGAAGT" .
node:102 a vg:Node ; vg:linksForwardToForward node:103 ; rdf:value "A" .
node:103 a vg:Node ; vg:linksForwardToForward node:104 ; rdf:value "TGA" .


alternative_node:103 a vg:Node ; vg:linksForwardToForward node:104 ; rdf:value "A" .
node:102 a vg:Node ; vg:linksForwardToForward alternative_node:103 .



```
In this example node:101 to node:104 are on the reference genome path.
They are also annotated to be a gene coding region.

```turtle
@prefix ex:<an  example ontology, could be SO> .

my_example:some_gene a ex:Gene, SO:Gene, vg:Path .
my_example:some_gene_step_1 a vg:Step ; vg:rank 1;vg:node node:101 .
my_example:some_gene_step_2 a vg:Step ; vg:rank 2;vg:node node:102 .
my_example:some_gene_step_3 a vg:Step ; vg:rank 3;vg:node node:103 .
my_example:some_gene_step_4 a vg:Step ; vg:rank 4;vg:node node:104 .

me:example:some_gene rdfs:seeAlso ENSEMBL:ESG00000XXXX . #and then pick up the annotation from the ENSEMBL gene build.
```

## Examples of using VG RDF

[2 ecoli genomes, with ensembl and uniprot annotation](/vgteam/vg/wiki/VG-RDF,-the-Ensembl-bacteria-E.-coli-genome-hack-attack)

## VG RDF limitations

At this moment VG RDF wants a fully embedded variation graph. e.g. all positions in vg json have a single edit which covers a whole node. This is to enable easy SPARQL queries where substring operations are rarely used.


## Annotations on pantograph and odgi

On top of VG RDF, we can describe the same path information on `odgi bin` and Pantograph format as well.

```turtle
# For representing odgi bin
bin:1 a odgi:Bin ; odgi:edge bin:2 ; rdf:value "ACT" .
bin:2 a odgi:Bin ; odgi:edge bin:3 : rdf:value "TGAAGT" .

path:1 a vg:Path .
path_step_1 a vg:Step ; odgi:bin bin:101 ; odgi:link bin:103 .
```