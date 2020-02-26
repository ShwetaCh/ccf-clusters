# ccf-clusters
Finding the number of cluster in a tumor using CCF values using GAP statistic (This statistic was created by Tibshirani to find optimal number of clusters for varied datasets https://statweb.stanford.edu/~gwalther/gap I tried to apply it on CCF to see how it works). This works only modestly with CCF for now, but can be optimized appropriately, I have done some tests to optimize the Kmax parameter which needs to be changed based on the platform. One can change the clustering method as well.

# Input
Maf with 'ccf_expected_copies_em' column

# Output
Number of clusters along with a CCF dot-plot, CCF densit-plot, and the Gap-statistic plot

