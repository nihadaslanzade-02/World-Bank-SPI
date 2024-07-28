# World-Bank-SPI
MATLAB analysis of the World Bank's Statistical Performance Indicators dataset, including PCA, K-means clustering, and various evaluation methods.

This repository contains MATLAB scripts and data files used for analyzing the Statistical Performance Indicators (SPI) for the year 2022. The analysis includes various data preprocessing, clustering, and dimensionality reduction techniques.

## Repository Contents

### MATLAB Scripts

#### analyze_spi_data.m

_Description:_ This script performs a detailed analysis of SPI data using different clustering techniques and selects the best performing K value based on explained variance.

_Usage:_ analyze_spi_data(X)

#### CDPCA.m

_Description:_ Performs clustering on both the rows (objects) and columns (variables) of the data matrix 𝑋, optimizing the explained variance of the data while ensuring the memberships matrices 𝑈 and 𝑉 are binary and row stochastic.
Vcdpca : Membership matrix for clustering variables
Ucdpca : Membership matrix for clustering objects
Acdpca : Loadings matrix (Principal Components)
Ycdpca : Projected data matrix
fcdpca : Final value of the objective function (explained variance)
incdpca : Number of iterations until convergence

_Usage:_ [Vcdpca,Ucdpca,Acdpca, Ycdpca,fcdpca,incdpca]=CDPCA(X, K, Q, Rndstart)

#### compute_explained_variance.m

Description: Computes the explained variance from PCA.
Usage: explained_variance = compute_explained_variance(data)

#### DKM.m

Description: Performs Double K-means clustering.
Usage: [centroids, assignments] = DKM(data, k)

#### elbow_method.m

Description: Uses the elbow method to determine the optimal number of clusters.
Usage: optimal_k = elbow_method(data, max_k)

#### FKM.m

Description: Performs Factorial K-means clustering.
Usage: [centroids, memberships] = FKM(data, k)

#### kmeansN.m

Description: Performs K-means clustering with a specified number of iterations.
Usage: [centroids, assignments] = kmeansN(data, k, n)

#### psF.m

Description: Computes Pseudo F-statistics.
Usage: scores = psF(data, n_components)

#### randPU.m

Description: Performs randomized projection and uplift.
Usage: proj_data = randPU(data, n_components)

#### REDKM.m

Description: Performs Reduced Dimensionality K-means clustering.
Usage: [centroids, assignments] = REDKM(data, k, n_components)
