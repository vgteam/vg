# VG in RDF

## Conceptual model

Nodes, Paths and Steps, are the three core parts of any VG graph in RDF.
Node in the VG RDF corersponds directly to the Node concept in the VG protobuf serialization.
Paths are a series of steps that represent a sequence of Node visits that are conceptually a 
biological sequence.
Steps are a way to easily find out the 



## Annotations on a VG

Path are annotations on the graph, some paths represent a reference genome, others are the genome of a patient.
Yet some represent a gene, regulatory element or any other feature that we can annotate on any biological sequence.

For example lets say we have a small VG graph:

```turtle
node:100 a vg:Node ; vg:linksForwardToForward node:101 ; rdf:value "ACT" .
node:101 a vg:Node ; vg:linksForwardToForward node:102 : rdf:value "TGAAGT" .
node:102 a vg:Node ; vg:linksForwardToForward node:103 ; rdf:value "A" .
node:103 a vg:Node ; vg:linksForwardToForward node:104 ; rdf:value "TGA"


alternative_node:103 a vg:Node ; vg:linksForwardToForward node:104 ; rdf:value "A"
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


## Linking a Genomic VG to a Protein VG


