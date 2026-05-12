# Hierarchical Statistical Analyzer

GUI-based Python framework for hierarchical and replicate-aware statistical analysis of biological datasets.

The software is designed for:
- stratified nonparametric testing
- blocked experimental designs
- hierarchical biological data
- replicate-aware inference
- skewed distributions commonly encountered in molecular biology

---

# Features

## Statistical methods

- Stratified Wilcoxon (van Elteren-style permutation test)
- Linear Mixed-Effects Models (LMM)
- Multiple testing correction:
  - Holm
  - Bonferroni
  - Benjamini-Hochberg FDR

## Experimental design support

- Biological replicate blocking
- Hierarchical data structures
- Planned contrasts
- Unequal sample sizes
- Skewed/non-Gaussian distributions

## GUI workflow

- CSV import
- Automatic column mapping
- Replicate assignment
- Contrast selection
- Exportable analysis results

---

# Installation

## Clone repository

```bash
git clone https://github.com/sebastienleonmichel/hierarchical-statistical-analyzer.git
cd hierarchical-statistical-analyzer
