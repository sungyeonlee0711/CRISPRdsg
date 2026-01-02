# CRISPRdsg
Prediction of sgRNA pair efficiency for CRISPR-mediated deletions.

## Description
This repository contains R scripts and input data used for:
- Feature correlation analysis
- Statistical comparison of sgRNA design features
- GLM and Random Forest modeling
- K-mer based modeling and feature importance
- Cross-validation, ROC/PR analysis, and Y-scrambling validation

## Files
- Analysis_pipeline.R: Full analysis pipeline
- Dual-sgRNA_design_input.txt: Design and basic variables input
- Kmer_countmatrix.txt: K-mer count matrix

## Requirements
- R (>= 4.0)
- Packages: dplyr, tidyr, randomForest, pROC, caret, PRROC, VennDiagram
