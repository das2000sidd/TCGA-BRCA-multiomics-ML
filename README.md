# TCGA-BRCA Multi-Omics Machine Learning

Deep-learning analysis of TCGA breast cancer data using omics data and clinical information to identify signature and candidate genes that can distinguish clinically relevant groups.

## Project overview

This project investigates whether transcriptomic and genomic features can be used to classify estrogen receptor (ER) status and identify molecular features associated with breast cancer outcomes.
Right now, the analysis only includes transcriptomic data

## Current results

- TCGA-BRCA RNA-seq and clinical data integrated for ML analysis
- PyTorch multilayer perceptron developed for ER-status classification
- ~0.85 macro ROC-AUC achieved on held-out test data
- Integrated Gradients and permutation importance used for model interpretation
- Candidate features evaluated using Cox proportional-hazards models

## Methods

- Python
- PyTorch
- pandas
- scikit-learn
- Captum / Integrated Gradients
- lifelines
- TCGA / GDC data

## Repository structure

- `Notebooks/` — analysis notebooks
- `Tables/` - important results
- `Figures/` — key figures
