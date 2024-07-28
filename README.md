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
_Description:_ Computes the explained variance for different K values in K-means clustering.
_Usage:_ explainedVariances = compute_explained_variance(X, K_values)

#### DKM.m
_Description:_ Performs Double K-means clustering on the data matrix.
_Usage:_ [Vdkm, Udkm, Ymdkm, fdkm, indkm] = DKM(X, K, Q, varargin)

#### elbow_method.m
_Description:_ Determines the optimal number of clusters using the pseudo F-statistic and the elbow method.
_Usage:_ elbow_method(X, K_values)

#### FKM.m
_Description:_ Performs Factorial K-means clustering on the data matrix.
_Usage:_ [Ufkm, Afkm, Yfkm, ffkm, infkm] = FKM(X, K, Q, Rndstart)

#### kmeansN.m
_Description:_ Performs K-means clustering with a specified number of clusters and random starts.
_Usage:_ [loopOtt, UOtt, fOtt, iterOtt] = kmeansN(X, K, Rndstart)

#### psF.m
_Description:_ Computes the pseudo F-statistic for clustering.
_Usage:_ [pf, Dw, Db] = psF(X, UOtt)

#### randPU.m
_Description:_ Generates a random partition of objects into classes.
_Usage:_ U = randPU(n, c)

#### REDKM.m
_Description:_ Performs Reduced Dimensionality K-means clustering on the data matrix.
_Usage:_ [Urkm, Arkm, Yrkm, frkm, inrkm] = REDKM(X, K, Q, Rndstart)


### Data and Plot Files

#### elbow_plot.png
_Description:_ A plot showing the elbow method for determining the optimal number of clusters.

#### explained_variances.png
_Description:_ A plot showing explained variance versus the number of clusters.

#### SPI_analysis.mat
_Description:_ Contains the results of the SPI analysis.

#### SPI_data1.xlsx
_Description:_ Contains the raw SPI data for 2022.


## Analysis Steps

### 1. Importing SPI Data
Import and load the SPI data into MATLAB for analysis.
### 2. Identifying Numeric and Categorical Columns
_Objective:_ Differentiate between numeric and categorical data for appropriate processing.
_Method:_ Use varfun, strcmp, and related functions to identify and separate numeric and categorical columns.
### 3. Converting Numeric Table to Array
_Objective:_ Convert the numeric data table into an array format suitable for mathematical operations.
_Method:_ Use table2array to perform the conversion.
### 4. Standardizing the Numeric Data
_Objective:_ Normalize the numeric data to have a mean of zero and a standard deviation of one.
_Method:_ Compute the mean (Xmean) and standard deviation (Xstd) of the data, and standardize the data using these values.
### 5. Computing Principal Component Analysis (PCA)
_Objective:_ Reduce the dimensionality of the data while retaining most of the variability.
_Method:_ Use the pca function to compute principal components.
### 6. Defining Initial K Value
_Objective:_ Establish an initial number of clusters for clustering methods.
_Method:_ Evaluate explained variances for different numbers of clusters (K) using PCA results.
### 7. Computing K-means Clustering
_Objective:_ Group the data into clusters based on the identified number of principal components and initial K value.
_Method:_ Use the kmeans function to perform clustering.
### 8. Determining Optimal K Value
_Objective:_ Identify the optimal number of clusters using both the elbow method and pseudo F-statistics.
_Methods:_
___Elbow Method:___ Use elbow_method.m to find the elbow point.
___Pseudo F-statistics:___ Use psF.m to compute the pseudo F-statistics, selecting the K with the highest pF.
### 9. Computing Reduced K-means
_Objective:_ Perform K-means clustering on reduced data based on initial and optimal K values.
_Method:_ Use REDKM.m to compute reduced K-means.
### 10. Computing Factorial K-means
_Objective:_ Execute factorial K-means clustering based on initial and optimal K values.
_Method:_ Use FKM.m for factorial K-means clustering.
### 11. Computing Clustering and Disjoint PCA
_Objective:_ Apply clustering and disjoint PCA based on initial and optimal K values.
_Method:_ Use CDPCA.m for the analysis.
### 12. Computing Double K-means
_Objective:_ Execute double K-means clustering based on initial and optimal K values.
_Method:_ Use DKM.m for double K-means clustering.
### 13. Analyzing and Comparing Results
_Objective:_ Compare the explained variances of different methods for different K values.
_Method:_ Use confusion matrices to evaluate the clustering performance across methods and K values.
### 14. Save and Display Analysis Results
_Description:_ Save the analysis results to a .mat file and display a summary of the key findings.
_Method:_ Use the save function to save results and disp to display the summary.

## Instructions for Use

Prerequisites
MATLAB R2022a or later
Statistics and Machine Learning Toolbox



Visualizing Results
Visualize the results using the provided plots elbow_plot.png and explained_variances.png.
Saving Results
Save the analysis results using save function in MATLAB.
save('spi_analysis_2022.mat', 'results', 'summary');


Contact

For any questions or issues, please contact [Your Name] at [Your Email].
