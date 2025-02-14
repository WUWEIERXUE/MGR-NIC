# MGR-NIC enables integrated cluster analysis of single-cell multi-omics data
Source codes and a demo of MGR-NIC are provided in this repository.
Running environment：``MATLAB R2019b`` or later. 
The external functions used in the MGR-NIC can be found in the ``MGR-NIC/`` folder.

###MGR-NIC includes the main functions below:
ClusteringMeasure_new.m:The ClusteringMeasure_new function calculates the clustering evaluation metrics to compare with the true_label.

DR_nmf.m:DR_nmf is a custom NMF function that returns the basis matrix W and the coefficient matrix H. The function is a custom NMF function that returns the base matrix W and the coefficients matrix H.

PMI.m：The PMI function calculates the Point Mutual Information (PMI) matrix M.

clustering_NIC.m：The clustering_NIC function performs multi-view clustering and returns the clustering result F and other variables.
clustering_mgrNIC.m：The clustering_mgrNIC function performs multi-view clustering and returns the clustering result F and other variables.

compute_f.m：Functions used to calculate precision, recall, and F1 scores.

constructW.m：The constructW function constructs the weight matrices W1, W2 for data X1 and data X2.

dataset1.mat：Simulation dataset 1.

dataset2.matSimulation dataset 2.

mESC.mat：A real single cell multiomics data utilized in the cell type clustering example. mESC for mouse embryonic stem cells (mESCs), including 13 cells cultured in 2i media and 64 serum-grown cells, which were profiled by parallel single-cell methylation and transcriptome sequencing technique scM&T-seq。

main_mgrNIC.m:Scripts using single-cell multi-omics data, showing how to run the code.

rand_index.m：Program for calculating the Adjusted Rand Index。
