# Test graph for content-aware node chopping

The graph is originally from the CHM13-based full HPRC v2.1 Minigraph–Cactus graph.

We extracted the largest snarl in chromosome 21 using gbz-base.
That snarl was between nodes 122616473 and 127206371.

Next we removed all but the first four paths from the graph and converted it to GBZ.
The result was the subgraph induced by those four paths.

Then we again extracted the largest snarl using gbz-base.
That snarl was between nodes 128885286 and 126090115.

The final graph is the unchopped version of that.
