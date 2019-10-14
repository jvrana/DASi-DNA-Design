# Algorithms

**Algorithm 1**: Given a DNA sequence and list of available fragments, primers,
and templates, determine a cost effective way to assembly given plasmid

1. create contigs of fragments, primers, and templates to the goal sequence
2. create assembly graph. Overhang. Golden-gate etc.
3. if cyclic, find shortest minimum cycles (modified floyd-warshall). If linear
find minimal shortest path between.

**Algorithm 2**: determine cost effective way to assembly sets of plasmids

