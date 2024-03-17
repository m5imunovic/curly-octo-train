The GFA file outputs edges using S label and L is an overlap of two edges (this is reverse to the OLC representation)
The string in GFA is kmer + truncated edge, so length of each string is > kmer. Coverage is integer coverage, in order
to get the true coverage this needs to be divided  by truncated edge length.

```cpp
size_t intCov() const { return cov; }
double getCoverage() const { return double(cov) / truncSize(); }
```

The GFA file only stores a single edge, its reverse complement needs to be induced from the contents of a file

The Dot file is more convenient, it already provides all edges with ids and this makes it easier to use together with
multiplicity info file. It stores edge as \[node_id_src -> node_id_dst\] and the additional edge id, truncated size and
normalized coverage as \[edge_id, size(coverage)\] label.
