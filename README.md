# TCGA-BRCA Multi-Omics Machine Learning

Deep-learning analysis of TCGA breast cancer data using omics data and clinical information to identify signature and candidate genes that can distinguish clinically relevant groups.

## Project overview

This project investigates whether transcriptomic and genomic features can be used to classify estrogen receptor (ER) status and identify molecular features associated with breast cancer outcomes.
Right now, the analysis only includes transcriptomic data

## Dataset

The dataset is uploaded onto zenodo here: https://zenodo.org/records/22953150?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjM1ZmFmNGM3LTM5ZjEtNGE5OC04YTY1LTk3ZmY1MGI2OGU0ZiIsImRhdGEiOnt9LCJyYW5kb20iOiJiM2NkMmI0MjE4Y2Y5ZmQ1ZjFmMDdmZGZkYzMxNThlYiJ9.LnrZUQ0E-xU1FWA9W1QYhJM_UeKQquf6CIJnBmQhUfemewYqeyd4w6-eyg8co9yh8KQ6pqk884nnnWwgKwA0eA

## Current results

- TCGA-BRCA RNA-seq and clinical data integrated for ML analysis
- PyTorch multilayer perceptron developed for ER-status classification
- ~0.85 macro ROC-AUC achieved on held-out test data
- Integrated Gradients and permutation importance used for model interpretation
- Candidate features evaluated using Cox proportional-hazards models
- The gene ENSG00000185008/ROBO2 which was among the genes with High permutation importance and -ve IG importance for class 0 (ER + ve) and +ve IG importance for class 1 (ER -ve) could be relevant to distinguish ER+ve from ER-ve breast cancers.  

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
