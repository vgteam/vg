# VG in RDF

## Conceptual model

`Node`s, `Path`s and `Step`s, are the three core parts of any VG graph in RDF.
A `Node` in the VG RDF corresponds directly to the Node concept in the VG protobuf serialization.
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

[2 ecoli genomes, with ensembl and uniprot annotation](https://github.com/vgteam/vg/wiki/VG-RDF,-the-Ensembl-bacteria-E.-coli-genome-hack-attack)

## VG RDF limitations

At this moment VG RDF wants a fully embedded variation graph. e.g. all positions in vg json have a single edit which covers a whole node. This is to enable easy SPARQL queries where substring operations are rarely used.


## Annotations on Pantograph

On top of VG RDF, we can describe the same path information on Pantograph format as well.

```ttl
<vg/zoom10> a vg:ZoomLevel ;
   vg:components <vg/zoom10/component1>, <vg/zoom10/component2> ;
   vg:zoomFactor 10 .
<vg/zoom1000> a vg:ZoomLevel ;
   vg:components <vg/zoom1000/component1>, <vg/zoom1000/component2> ;
   vg:zoomFactor 1000 . // zoomFactor is binWidth here.
<vg/zoom1000/component1> a vg:Component ;
   vg:componentRank 1 ;   # The order of component is inferred by rank.
   vg:forwardComponentEdge <vg/zoom1000/component2> ;
   vg:bins <vg/zoom1000/component1/bin1>, <vg/zoom1000/component1/bin2> .
<vg/zoom1000/component1/bin1> a vg:Bin ;
   vg:forwardBinEdge <vg/zoom1000/component2/bin2> ;
   vg:binRank 1 ;
   vg:cells <vg/zoom1000/component1/bin1/cell1>, <vg/zoom1000/component1/bin1/cell2> .
<vg/zoom1000/component2/bin2> a vg:Bin ;
   vg:reverseBinEdge <vg/zoom1000/component2/bin3> ;
   vg:forwardBinEdge <vg/zoom1000/component2/bin3> ;
   vg:binRank 2 ;
   vg:cells <vg/zoom1000/component1/bin2/cell1>, <vg/zoom1000/component1/bin2/cell2> .
<vg/zoom1000/component1/bin1/cell1> a vg:Cell ;
   vg:positionPercent 0.04 ;
   vg:inversionPercent 0.98 ;
   vg:cellRegion <path1/region/6-100> .  # To infer firstNucleotide and last Nucleotide. faldo:begin of stepRegion is the first position. faldo:end of cellRegion is the last position.
<vg/zoom1000/link1> a vg:Link ; # This is a non-linear connection between Bins.
   vg:arrival <vg/zoom1000/component1/bin1> ;
   vg:departure <vg/zoom1000/component2/bin2> ;
   vg:forwardLinkEdge <vg/zoom1000/link2> ;
   vg:linkRank 1 ;
   vg:linkPaths <path1> <path2> ; # Participants of the link
   vg:linkZoomLevel <vg/zoom10> .
<path1/region/6-100> a faldo:Region ;
   faldo:begin <path1/position/6>  ;
   faldo:end <path1/position/100>  .
```
