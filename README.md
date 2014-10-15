# Tree length project analyses

## Sebastian Duchene 15 Oct

The assignment of trees to topology clusters is in clust_k_10.txt. Use this file to concatenate the alignments and estimate an ML tree topology. Then optimise the branch lengths for each gene. This should result in 10 files which would contain the gene trees. Run clockstar on each of these files. This will give the pacemakers for each number of clusters.

To determine the ranking of tree according to their length, load the ML trees in R and add the branch lengths, then save this to a file... 

The script tree_names.py can be executed from python or terminal to find whether there are duplicate trees in a .trees file (newick with tree names).
