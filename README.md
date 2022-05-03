# mds_latent

- This repository stores the code required to cluster binary mutation profiles of MDS cases using a single-layer PCA-like autoencoder framework. 'pca_cluster.py' requires a binary matrix with rows as observations and columns as features and generates a consensus-matrix 'pca_cc_res.npy' quantifies the frequency of pairwise-clusterings.
- Additional scripts are used for preprocessing and plotting; preparing the input data (process_mds.R), imputation of missing data (rf_impute.py), generating figures (generate_plots.R)
- Checkout https://drmz.shinyapps.io/mds_latent/ for a simple implementation of predictive tool for classification
