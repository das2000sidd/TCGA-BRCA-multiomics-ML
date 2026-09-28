# TCGA-BRCA Multi-Omics Machine Learning

Deep-learning analysis of TCGA breast cancer data integrating transcriptomic, mutational and copy-number features to investigate molecular predictors of estrogen receptor (ER) status.

## Project overview

This project evaluates whether multi-omics data from **TCGA-BRCA** can be used to classify estrogen receptor (ER) status using neural-network models.

The analysis progressively evaluates:

* RNA-seq gene expression
* Somatic mutation features
* Copy-number variation (CNV)
* RNA + mutation integration
* RNA + CNV integration
* Model interpretation using permutation importance and Integrated Gradients
* Association of candidate genes with clinical outcomes using survival analysis

The primary prediction task is binary ER-status classification:

* **Class 0: ER-positive**
* **Class 1: ER-negative**

The project is designed to evaluate both **predictive performance** and whether additional omics modalities provide complementary information beyond transcriptomic data.

---

## Dataset

The analysis uses TCGA-BRCA data obtained through the Genomic Data Commons (GDC).

Processed datasets used in the analysis are available through Zenodo:

[Zenodo dataset](https://zenodo.org/records/22953150?utm_source=chatgpt.com)

### Dataset summary

| Data type                |                   Features / samples |
| ------------------------ | -----------------------------------: |
| RNA-seq                  |            919 samples × 3,796 genes |
| Somatic mutation         |                394 mutation features |
| CNV                      | 1,093 samples × 3,765 genes after QC |
| RNA + CNV common cohort  |                          917 samples |
| Independent RNA test set |                          276 samples |

---

## Analysis workflow

```text
                    TCGA-BRCA
                        │
        ┌───────────────┼────────────────┐
        │               │                │
      RNA-seq        Mutation           CNV
        │               │                │
   3,796 genes      394 features     Gene-level CNV
        │               │                │
        │        Univariate AUC      Correlation
        │        + MLP analysis      clustering
        │               │                │
     OOF permutation importance          │
        │                                │
     Top 500 genes                  CNV clusters
        │                                │
        └───────────────┬────────────────┘
                        │
                  MLP classification
                        │
                    ER status
                        │
              Model interpretation
                        │
             Survival analysis
```

---

# Results

## 1. RNA expression

A multilayer perceptron (MLP) was trained using RNA-seq expression data.

Feature selection was performed using **5-fold out-of-fold permutation importance** within the development cohort. Models using 100, 200 and 500 genes were subsequently evaluated on the validation set.

| RNA model        | Features | Train n | Validation n |   ROC-AUC |  Accuracy |
| ---------------- | -------: | ------: | -----------: | --------: | --------: |
| Full RNA         |    3,796 |     514 |          129 |     0.963 |     0.938 |
| Selected RNA     |      100 |     514 |          129 |     0.961 |     0.922 |
| Selected RNA     |      200 |     514 |          129 |     0.963 |     0.907 |
| **Selected RNA** |  **500** | **514** |      **129** | **0.971** | **0.938** |

The 500-gene representation was selected for subsequent multimodal experiments.

### Independent test-set performance

After feature selection and model development, the final 500-gene RNA model was retrained on the development data and evaluated on an **independent held-out test set (n = 276)**.

| Metric      | Performance |      95% CI |
| ----------- | ----------: | ----------: |
| ROC-AUC     |   **0.942** | 0.896–0.979 |
| Accuracy    |   **94.6%** |  91.7–97.1% |
| Sensitivity |   **88.9%** |  80.3–96.1% |
| Specificity |   **96.2%** |  93.6–98.6% |
| PPV         |   **87.5%** |  78.8–95.1% |
| NPV         |   **96.7%** |  94.0–99.0% |

Confusion matrix:

```text
                 Predicted
                 ER+   ER-
Actual ER+       205     8
Actual ER-         7    56
```

The independent test set was not used for feature selection, model selection or preprocessing parameter estimation.

---

## 2. Mutation features

A mutation-only MLP showed more limited predictive performance:

| Model         | Features | Validation ROC-AUC | Accuracy |
| ------------- | -------: | -----------------: | -------: |
| Mutation-only |      394 |              0.663 |    0.775 |

However, individual mutation features contained substantial biological signal.

### TP53

TP53 mutation showed the strongest univariate association with ER status among the evaluated mutation features.

| Feature  | Direction-independent AUC |
| -------- | ------------------------: |
| **TP53** |                 **0.809** |
| PIK3CA   |                     0.624 |
| GATA3    |                     0.584 |
| CDH1     |                     0.573 |
| AR       |                     0.559 |
| MAP3K1   |                     0.556 |

TP53 mutation frequency:

* ER-positive: **20.5%**
* ER-negative: **82.2%**

This illustrates that individual genomic features can be highly informative even when the mutation-only neural-network model has relatively modest overall predictive performance.

Adding TP53 to the 500-gene RNA model gave:

| Model      | Validation ROC-AUC | Accuracy |
| ---------- | -----------------: | -------: |
| RNA only   |              0.971 |    0.938 |
| RNA + TP53 |              0.968 |    0.938 |

---

## 3. Copy-number variation

CNV features were evaluated independently before integration with RNA expression.

The full gene-level CNV model achieved:

**ROC-AUC = 0.942**

Individual-gene permutation selection substantially reduced performance:

| CNV representation | Validation ROC-AUC | Accuracy |       |
| ------------------ | -----------------: | -------: | ----- |
| Full CNV           |        3,764 genes |    0.942 | 0.922 |
| Top 100 genes      |                100 |    0.880 | 0.891 |
| Top 200 genes      |                200 |    0.856 | 0.853 |
| Top 500 genes      |                500 |    0.879 | 0.868 |

This suggested that predictive CNV information was distributed across correlated genomic features rather than concentrated in a small number of individual genes.

### Correlation-based CNV representation

Genes were clustered using correlation-based distances based on:

```text
distance = 1 - |correlation|
```

Cluster-level CNV features were then generated by averaging CNV values across genes within each cluster.

This representation improved validation performance to:

**ROC-AUC = 0.955**
**Accuracy = 0.938**

---

## 4. RNA + CNV integration

RNA and CNV were evaluated using the same 917-sample cohort and identical train/validation/test splits.

A matched RNA-only model was trained using the same 512 training samples and 129 validation samples.

| Model        | Features               | Train n | Validation n |   ROC-AUC |  Accuracy |
| ------------ | ---------------------- | ------: | -----------: | --------: | --------: |
| **RNA only** | 500 RNA genes          |     512 |          129 | **0.981** | **0.961** |
| RNA + CNV    | 500 RNA + CNV clusters |     512 |          129 |     0.960 |     0.938 |

On this validation cohort:

* RNA-only model: **124/129 correctly classified**
* RNA + CNV model: **121/129 correctly classified**
* RNA-only correct → RNA+CNV incorrect: **5 samples**
* RNA-only incorrect → RNA+CNV correct: **2 samples**

Thus, although CNV showed substantial predictive signal independently, adding CNV through the current feature-fusion approach did **not improve the RNA model**.

This does not establish that CNV contains no complementary biological information; rather, the current CNV representation and fusion architecture did not provide additional predictive benefit in this validation experiment.

---

# Overall model comparison

| Modality / model     |                   Features | Train n | Validation n |   ROC-AUC |  Accuracy |
| -------------------- | -------------------------: | ------: | -----------: | --------: | --------: |
| RNA                  |                3,796 genes |     514 |          129 |     0.963 |     0.938 |
| RNA                  |              Top 100 genes |     514 |          129 |     0.961 |     0.922 |
| RNA                  |              Top 200 genes |     514 |          129 |     0.963 |     0.907 |
| **RNA**              |          **Top 500 genes** | **514** |      **129** | **0.971** | **0.938** |
| Mutation             |               394 features |     514 |          129 |     0.663 |     0.775 |
| RNA + mutation       |     500 RNA + 394 mutation |     514 |          129 |     0.963 |     0.938 |
| RNA + TP53           |             500 RNA + TP53 |     514 |          129 |     0.968 |     0.938 |
| CNV                  |                3,764 genes |     512 |          129 |     0.942 |     0.922 |
| CNV                  |              Top 100 genes |     512 |          129 |     0.880 |     0.891 |
| CNV                  |              Top 200 genes |     512 |          129 |     0.856 |     0.853 |
| CNV                  |              Top 500 genes |     512 |          129 |     0.879 |     0.868 |
| **CNV clusters**     |   **Correlation clusters** | **512** |      **129** | **0.955** | **0.938** |
| **RNA + CNV**        | **500 RNA + CNV clusters** | **512** |      **129** | **0.960** | **0.938** |
| **Matched RNA-only** |                **500 RNA** | **512** |      **129** | **0.981** | **0.961** |

> Validation results are reported for model development and modality comparison. The independent 276-sample test set was reserved for the locked RNA model.

---

# Model interpretation

Two complementary approaches are being used to identify features contributing to model predictions:

### Permutation importance

Feature importance was estimated using out-of-fold permutation importance to reduce dependence on a single train/validation split.

### Integrated Gradients

Integrated Gradients was used to estimate feature contributions to the neural-network output.

Several genes showed evidence of importance across the two approaches. One example is **ROBO2 (ENSG00000185008)**, which showed high permutation importance and class-specific Integrated Gradients attribution.

These features are considered **candidate predictive features**, rather than causal biomarkers.

---

# Survival analysis

Candidate genes identified from the classification models are additionally evaluated using **Cox proportional-hazards models** to investigate whether they are associated with clinical outcomes.

Classification performance and survival association are treated as separate analyses: a gene can contribute to ER-status prediction without necessarily being associated with survival.

---

# Key findings

1. **RNA expression provided strong predictive performance for ER-status classification.**
2. A 500-gene RNA representation retained strong performance while substantially reducing dimensionality from 3,796 genes.
3. The final locked RNA model achieved **ROC-AUC 0.942 (95% CI 0.896–0.979)** on an independent test set.
4. Mutation data showed weaker overall predictive performance, although individual mutations such as **TP53** were strongly associated with ER status.
5. CNV showed substantial standalone predictive signal, with correlation-based CNV clusters achieving **validation ROC-AUC 0.955**.
6. Adding CNV to the RNA model did **not improve performance** in the matched validation experiment.
7. Model interpretation identified candidate genes that can be investigated further using biological annotation and survival analysis.

---

# Methods

* Python
* PyTorch
* pandas
* NumPy
* scikit-learn
* SciPy
* Captum / Integrated Gradients
* lifelines
* TCGA / GDC data

### Machine-learning approach

* Multilayer perceptron (MLP)
* Stratified train/validation/test splitting
* Standardization using training data only
* Out-of-fold permutation feature importance
* Correlation-based CNV feature clustering
* Independent held-out test evaluation
* Bootstrap confidence intervals

---

## Project status

**Current stage:** RNA, mutation and CNV analyses completed. DNA methylation analysis is the next planned modality for evaluation.

The project is ongoing, with future work focused on evaluating whether additional molecular modalities provide complementary predictive information and on biological interpretation of model-derived candidate features.
