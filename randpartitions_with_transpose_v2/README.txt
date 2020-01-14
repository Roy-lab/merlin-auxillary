Usage: 
./makePartitions inputdata partitions outputdir partitionsize partitiontype[rand|exclusive] orientation[normal|transpose]

./makePartitions example_in/geneexp.txt 100 example_out/randpartitions_n20/ 20 rand normal

Input expression file is in example_in/geneexp.txt
	Should not have sample names (just genes and expression)
partitions: is the number of partitions to make
outputdir: is the output directory to store all files (this directory needs to be created before running this program)
partitionsize: size of the partition, that is number of samples
partitiontype: rand will do random subsampling without replacement and will produce datasets that are potentially overlapping with each other. exclusive will generate non-overlapping partitions.
orientation: normal means genes are on rows, samples on columns. transpose means genes are on columns, samples are on rows. Input and output will be in the same orientation.


UPDATE 3/26/19 - Viswesh Periyasamy
created partitionRandWithSamples.py which accepts gene expression and writes subsampled files with sample labels
NOTE: does not accept exclusive as an option, just does random sampling without replacement
NOTE: requires pandas and numpy environment

usage:
python3 partitionRandWithSamples.py inputdata partitions outputdir partitionsize orientation[normal|transpose]

python3 partitionRandWithSamples.py test_with_samples.txt 10 example 3 normal

UPDATE 11/27/19 - Matt Stone
Added makePartitions.py with more control over outputs. 
Other added features:
    - can auto-compute partition size as fraction of total samples
    - can read gzipped input
    - optionally gzip output
    - specify random seed when sampling for reproducibility
