# Hierarchical Statistical Analyzer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19500564.svg)](https://doi.org/10.5281/zenodo.19500564)

GUI-based Python framework for hierarchical and replicate-aware statistical analysis of biological datasets.

**Current version: 4.0.0.** This version changes the linear mixed-effects model (LMM) to include fixed condition and categorical experiment effects, with a random intercept for each experiment × condition sample.

## What's new in v4.0.0

- Revised LMM structure for repeated observations within biological samples.
- Explicit B-minus-A condition coding and extraction of the condition p-value by name.
- Checks for insufficient replication, lack of repeated observations, confounded condition/experiment effects, nonconvergence, and invalid condition p-values.
- LMM appears first in the analysis options and exported results, with its model structure displayed in the GUI and result method label.

**Migration from v3:** v3 used `response ~ group` with a random intercept for `replicate_id`. The v4 model changes the statistical assumptions, so LMM p-values can differ for the same dataset. Review the mapping below before rerunning an analysis. The permutation tests, response-transform choices, Holm correction, and CSV column names are retained.

## Statistical methods

- Linear mixed-effects models fitted separately for each planned pairwise contrast.
- Global stratified rank permutation test across groups included in the selected contrasts.
- Pairwise stratified Wilcoxon (van Elteren-style) permutation tests.
- Within-experiment permutations that preserve group sizes; the GUI uses 10,000 permutations.
- Holm correction across valid planned pairwise p-values, separately within each statistical engine.
- Block-aware Cliff's delta, group medians, and median differences.

### LMM model and data mapping

For a contrast **A vs B**, the implemented model is:

```text
response ~ condition_B + C(replicate_id)
random intercept: experiment × condition
```

`condition_B` is 0 for A and 1 for B. `replicate_id` is treated as a categorical experiment identifier. All observations with the same experiment and condition share one random intercept.

| Data level | Meaning in v4 |
| --- | --- |
| Observation | An individual measurement, such as one cell or DNA fiber |
| Group / condition | The treatment or biological group being compared |
| Replicate ID | The experiment/block shared by conditions measured in the same experiment |
| Experiment × condition | One biological sample containing repeated observations |

Use the same Replicate ID for paired conditions from the same experiment. For example, map control and treatment from experiment 1 to `exp1`, and both conditions from experiment 2 to `exp2`.

The model assumes **one biological sample per experiment × condition combination**. Multiple columns mapped to that same combination are pooled under the same sample intercept. Independent cultures within the same combination cannot be represented as separate biological samples by this mapping.

The LMM requires repeated observations within samples and at least two experiments containing both conditions for each contrast. It rejects designs where condition and experiment effects cannot be separated. These are minimum fit requirements, not a guarantee of adequate replication or reliable inference.

Fits use maximum likelihood (`reml=False`) and two-sided asymptotic Wald p-values for the condition coefficient. There is no small-sample correction. Failed fits remain visible as `LMM failed: ...` rows with missing p-values and blank significance symbols. The script checks convergence and p-value validity; it does not automatically reject every singular-covariance or boundary-variance warning.

### Response transforms

Available LMM transforms are **None, log2, ln, and log10**. They affect only the LMM response. Permutation tests and descriptive effect summaries use the original measurement scale.

Log transforms require strictly positive responses. Nonpositive values produce a failed LMM result for the affected contrast; they are not silently discarded or replaced with a pseudocount. Zero and negative values remain available for untransformed analyses. Blank, nonnumeric, and nonfinite input entries are excluded during reshaping.

## Installation

Clone the repository and install its dependencies:

```bash
git clone https://github.com/sebastienleonmichel/hierarchical-statistical-analyzer.git
cd hierarchical-statistical-analyzer
python -m pip install -r requirements.txt
```

Dependencies are NumPy, pandas, SciPy, statsmodels, and openpyxl. No additional dependency is introduced in v4.0.0.

Tkinter is included with most standard Python installations. On systems where it is packaged separately, install the appropriate Tk package for your Python distribution.

## Usage

```bash
python hierarchical_statistical_analyzer.py
```

1. Open a wide-format `.csv` or `.xlsx` file, with measurements arranged down each input column.
2. Review the automatically suggested group and Replicate ID mapping for every column.
3. Select the planned pairwise contrasts.
4. Select the LMM, stratified Wilcoxon, or both, and choose the LMM response transform.
5. Run the analysis. Background execution reports status and elapsed time; results appear on screen and are saved beside the input file as `<input-name>_statistical-analysis.csv`.

## Results

When both engines are selected, LMM rows come first, followed by the global and pairwise permutation results.

The exported columns are:

```text
Contrast, Method, Response_transform, p_raw, p_holm,
Cliffs_delta_B_minus_A, Median_A, Median_B,
Median_difference_B_minus_A, Symbol
```

Positive effect summaries mean B tends to be higher than A. Cliff's delta compares observations within complete blocks and weights blocks by their number of A–B observation pairs. `Median_A` and `Median_B` are medians of block-level medians; the median difference is the median of within-block B-minus-A differences. These are descriptive summaries on the original scale, not fitted LMM coefficients.

`Symbol` is the final column and is based only on `p_holm`: `****` ≤ 0.0001, `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05, and `ns` otherwise. Missing adjusted p-values give a blank symbol. The global test reports its raw p-value and is excluded from pairwise Holm correction.

## Citation

Use the version information in [`CITATION.cff`](CITATION.cff) and identify the release used in your analysis. The DOI badge links to the Zenodo record for **all versions**; for a reproducible citation, select the DOI for your exact archived version when available.

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE).
