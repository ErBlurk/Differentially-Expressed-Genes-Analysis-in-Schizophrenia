# Differentially Expressed Genes Analysis in Schizophrenia

## Overview

This repository contains code and resources for La Sapienza bioinformatics project (course 1047220, A.A. 2025), aimed at identifying differentially expressed genes (DEGs) associated with schizophrenia. The analysis utilizes publicly available gene expression data from the Gene Omnibus to uncover potential biomarkers and understand the molecular underpinnings of the disorder.

## Data Source

The analysis is based on gene expression data obtained from the Gene Expression Omnibus (GEO), specifically the [GSE38484](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38484) dataset. 

## Project Structure

- **BioMart/**: Scripts for retrieving gene annotations and related information using the BioMart database.
- **PCA/**: Principal Component Analysis scripts to assess sample clustering and variability.
- **Results/**: Output files including DEG lists, plots, and enrichment analysis results.
- **enrichR/**: Scripts for performing gene set enrichment analysis using the Enrichr tool.
- **DEG_AnalysisGEO.R**: Script to perform differential expression analysis on the downloaded data.
- **GEOquery.R**: Main script for downloading and preprocessing data from the GEO database.
- **NCBIquery.R**: Scripts for querying NCBI databases for additional gene information.

## Methodology

This project applies a bioinformatics workflow to identify differentially expressed genes (DEGs) in schizophrenia using whole-blood data from the GSE38484 dataset.

1. Gene expression data were imported from GEO and annotated with clinical information (control vs. schizophrenia).
2. Low-expression and low-variability genes were removed after averaging, log-transformation, and filtering by interquartile range (IQR).
3. Differential expression was calculated using log fold-change and t-tests, followed by p-value adjustment with the False Discovery Rate (FDR).
4. Genes with low fold-change or non-significant adjusted p-values were excluded, and volcano plots were used to visualize the results.
5. A heatmap and PCA plot were generated to explore clustering patterns and expression differences between groups.
6. Significant DEGs were exported to a text file including log fold-change and adjusted p-values.
7. Functional enrichment analysis was performed using Enrichr to identify biological processes and pathways affected in schizophrenia.
8. Additional analysis was done using miRTarBase and Mienturnet to identify key regulatory miRNAs
9. Finally a literature-based research has been performed to review the biological relevance of the results.

## Outputs

Processed data and visualizations are stored in the `Results/` directory. Highlights include:

- Lists of upregulated and downregulated genes in schizophrenia.
- PCA plots illustrating sample clustering.
- Enrichment analysis results pointing to dysregulated pathways.

