# merlin-downstream
Downstream analysis for MERLIN inferred networks

# Making sub-samples
`randpartitions_with_transpose` is used to create sub-samples.
You can run it on an expression matrix with:
* genes on rows and samples on columns (normal) or 
* samples on rows and genes on columns (transpose).
It assumes there is no header file. 
So:
* For normal, first column is gene names, there is no sample name in first line
* For transpose, first row is gene names, there is no sample name in first column
For example:
```
./makePartitions trans.txt 100 outdir/ 50 rand transpose
```
* trans.txt is a transposed expression matrix
* it creates 100 subsamples (dataset0.txt to dataset99.txt)
* in outdir/ 
* each with 50 randomly selected samples.

# Making consensus networks
`estimateedgeconf` creates a consensus network from a list of network files:
```
./estimateEdgeConf network_files.txt 0 output_net_ alledges
```
where `network_files.txt` has the location of individual network files and the output (`output_net_alledge.txt`) will contain edges, and percentage of times the edges were seen in individual networks (1 means 100%).

# Making consensus modules
`assessclusterstab` creates a co-clustering matrix that shows how many times two genes were in the same module.
 
```
./assessClusterStab module_files.txt sims.txt
```
where `module_files.txt` has the list of module assignments files (one per line) and sims.txt will be the co-clustering matrix.

`optimalleaforder` applies hierarchical clustering to co-clustering matrix from the previous step.

```
./reorder sims.txt matrix consensus_module_0.3 0.3
```
where 0.3 is the threshold for creating the hierarchical modules.

