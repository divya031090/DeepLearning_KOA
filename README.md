#Introduction

This repository contains code for a comprehensive analysis framework utilizing multiomic multimodal data, including clustering using Variational Autoencoders, differential expression analysis, and integrative machine learning for predicting longitudinal response.

##Repository Structure

Clustering_Variational_Autoencoder.py: This Python script implements a Variational Autoencoder for clustering data based on integrating multimodal data.

Diff_Expr_signature.R: An R script for conducting differential expression analysis, identifying significant differences between groups.

Std_mean_diff_signature.R: An R script for calculating standardized mean differences, used for signature-based analysis.

Signature_plots.R: This R script generates signature plots for the top 10 variables identified wihtin each cluster.

Diff_expr_response.R: An R script for selecting variables for machine learning analysis based on differential expression w.r.t longitudinal response.

Integrative_ML.R: An R script that outlines an integrative machine learning framework to predict longitudinal response using the multimodal data domains.

##Getting Started

##Prerequisites
Python 3.10
R 4.2
Required Python libraries: numpy, pandas, tensorflow, keras, scikit-learn
Required R packages: DESeq2, limma, ggplot2, dplyr, caret

##Usage
Run the files in this order:
Clustering with Variational Autoencoder
Run the Clustering_Variational_Autoencoder.py script to perform clustering on your dataset.

python Clustering_Variational_Autoencoder.py

Differential Expression Analysis
Use the Diff_Expr.R script to perform differential expression analysis for the signature within each cluster

Rscript Diff_Expr.R

Standardized Mean Differences
Calculate standardized mean differences using std_mean_diff.R for the signature within each cluster:

Rscript std_mean_diff.R

Signature Plots
Generate signature plots for the top 10 variables identified as significant through both Standardized Mean Differences and Differential Expression Analysis

Rscript signature_plots.R

Variable Selection for ML Analysis
Pre-process and select variables for machine learning modeling based on abs(log fold change) of 1.1 for metabolites and 1.5 for miRNAs:

Rscript Diff_expr_response.R

Integrative Machine Learning Framework
Predict longitudinal response using the integrative machine learning framework 

Rscript Integrative_ML.R

License

This project is licensed under the MIT License - see the LICENSE file for details.
