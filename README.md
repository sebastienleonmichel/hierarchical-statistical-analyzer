# Hierarchical Statistical Analyzer — ZetaStats

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19500564.svg)](https://doi.org/10.5281/zenodo.19500564)

GUI and command-line tool for hierarchical, replicate-aware analysis of biological measurements. **Current version: 5.1.0**, displayed as **ZetaStats v5.1.0** in the application.

Python handles input, mapping, the interface, and Excel export. R's **lme4**, **lmerTest**, and **jsonlite** provide one global linear mixed model, fitted by REML, with Satterthwaite degrees of freedom for Type III F tests and planned contrast t tests.

## What's new in v5.1.0

- Explicit **Full**, **Reduced**, and **Custom** random-effects choices in the GUI and Python API; full remains the default.
- CLI presets with `--random-effects full|reduced` and individual `--drop-replicate-condition` / `--drop-replicate-treatment` switches.
- Experiment and biological sample random intercepts remain in every selected structure; fixed effects and the global-fit workflow are unchanged.
- Exports record the selected structure, exact formula, included/excluded/not-applicable terms, and model status, including when inference is withheld.
- Scrollable analysis options and output filenames that identify reduced or custom analyses.

Choose the random-effects structure for the experimental design before interpreting results. Changing structure can change estimates, degrees of freedom, and p-values; v5.1.0 does not select a simpler model automatically after a singular or failed fit.

### Features introduced in v5.0.0

- Separate condition and treatment factors, including their fixed interaction.
- One global model fitted to all included samples; selected contrasts use that same fit.
- Random effects for experiment, experiment × condition, experiment × treatment, and biological sample.
- Optional normalization to the mean of reference replicate means, followed by an optional response transform.
- Three-sheet XLSX export with descriptive statistics, variance components and Type III tests, and Holm-adjusted contrasts.
- Inference withheld for singular, nonconverged, or otherwise invalid fits; no automatic simpler-model fallback.
- Command-line analysis and backend checks in addition to the GUI.

**Migration from v4:** this is a change in statistical model, input mapping, dependencies, and output format. v4 fitted a separate Python/statsmodels LMM for each contrast and also offered stratified permutation tests. v5 uses the R-backed global model and does not include those permutation tests or Cliff's delta summaries. Review the sample mapping before rerunning existing data; v4 and v5 p-values are not interchangeable. Earlier versions remain available in the repository's release history.

## Installation

Use Python 3.10 or newer with Tkinter, and install the Python dependencies:

```bash
git clone https://github.com/sebastienleonmichel/hierarchical-statistical-analyzer.git
cd hierarchical-statistical-analyzer
python -m pip install -r requirements.txt
```

The Python dependencies are NumPy, pandas, and openpyxl. SciPy and statsmodels are no longer used. Tkinter is included with many Python installations; install the Tk package for your Python distribution if needed.

Install [R](https://cran.r-project.org/) separately, then install the backend packages from an R console:

```r
install.packages(c("lme4", "lmerTest", "jsonlite"), repos = "https://cloud.r-project.org")
```

Check the backend:

```bash
python hierarchical_statistical_analyzer.py --check-backend
```

The script searches for `Rscript` on PATH and at common macOS installation locations. Set `ZETASTATS_RSCRIPT` to an executable path or use `--rscript /path/to/Rscript` when needed. These source files require a separate R installation; a self-contained application bundle is not included here.

## Input and biological units

Load a wide-format `.csv` or `.xlsx` file with **one biological sample per column** and individual measurements down the rows. The GUI reads the first Excel worksheet; the command line can select a worksheet by name with `--sheet`.

Suggested column names follow `Condition_treatment_replicateID`, for example:

```text
WT_0J_N1, WT_2J_N1, SCAIko_0J_N1, SCAIko_2J_N1,
WT_0J_N2, WT_2J_N2, SCAIko_0J_N2, SCAIko_2J_N2
```

Names are parsed from the right, so condition names may contain underscores. The GUI allows all three labels to be edited and annotation columns to be excluded.

| Data level | Meaning |
| --- | --- |
| Observation | One measurement, such as a cell or DNA fiber |
| Condition | A categorical biological factor, such as genotype |
| Treatment | A categorical treatment factor, such as UV dose |
| Replicate ID | The independent experiment shared by its matched samples |
| Condition × treatment × replicate | One biological sample with repeated observations |

Use the same replicate ID for samples measured within the same independent experiment. Cells remain nested within samples; their count does not increase biological replicate n. Multiple columns cannot map to the same condition–treatment–replicate combination. Combine technical measurements explicitly before importing. The mapping does not represent multiple independent cultures within one such combination.

Blank cells are omitted. Duplicate headers, nonnumeric or infinite entries in included columns, empty included samples, and duplicate sample mappings are rejected.

## Statistical model

For two factors with more than one level, the default **Full** model is:

```r
response ~ condition * treatment +
  (1 | replicate_id) +
  (1 | replicate_id:condition) +
  (1 | replicate_id:treatment) +
  (1 | sample_id)
```

`sample_id` identifies condition × treatment × replicate. Condition, treatment, and their interaction are fixed effects. When included, random replicate × treatment deviations allow treatment responses to vary between experiments using a shared categorical variance component; this is not an unrestricted correlated random-slope model.

The **Reduced** model explicitly excludes both optional replicate × factor terms:

```r
response ~ condition * treatment + (1 | replicate_id) + (1 | sample_id)
```

A **Custom** structure includes only one of `(1 | replicate_id:condition)` and `(1 | replicate_id:treatment)`. Experiment and biological sample intercepts are mandatory in all structures. All tests use the selected global fit. The selected structure and exact formula are recorded in `2_Variance_TypeIII`, together with each optional term's requested and applied state. Aliasing checks apply to the selected random components.

If one factor has a single level, the fixed part contains the other factor and the random part contains replicate and sample intercepts. The optional replicate × factor components are then marked not applicable, regardless of their requested setting. Treatment is categorical, including when its labels are numeric doses; the model does not fit a continuous dose-response slope.

The design must include at least two condition–treatment combinations, every combination of the included factor levels, at least two independent replicates overall and per combination, and at least two observations per sample. Aliased random components and rank-deficient fixed effects are rejected. These are minimum fit requirements, not a guarantee of adequate replication. With fewer than five experiments, the exported metadata flags that variance estimates and denominator degrees of freedom may be imprecise.

Fits use REML with sum contrasts. Type III fixed-effect F tests and planned B-minus-A contrast t tests use Satterthwaite degrees of freedom from **lmerTest**. All included observations inform the same fit, even when only a subset of contrasts is selected.

Holm adjustment applies across **all selected contrasts**, including failed tests in the planned family size. Type III omnibus tests are reported separately with unadjusted p-values. The 95% contrast confidence intervals are individual intervals and are not adjusted for multiplicity.

Singular or boundary variance fits, convergence failures, invalid covariance, and invalid individual tests withhold the corresponding inferential results. Descriptive output and available fit details remain visible with an explicit status. No simpler model is substituted automatically.

### Normalization and transforms

Available response transforms are **None, log2, ln, and log10**. Logs require strictly positive values; no pseudocount is added.

Optional normalization divides every observation by one shared positive divisor: the arithmetic mean of the reference condition–treatment group's replicate means, with equal weight for each reference replicate. Normalization occurs **before** transformation. Inference is conditional on this observed divisor; its uncertainty is not propagated as an additional model term.

The fitted model assumes Gaussian residuals with a common variance on the response scale. A log transform may help but does not guarantee these assumptions. Review the exported conditional residual summaries and variance patterns.

## Usage

Open the GUI:

```bash
python hierarchical_statistical_analyzer.py
```

1. Open a CSV or XLSX file.
2. Review condition, treatment, and replicate mapping; exclude annotation columns.
3. Click **Update contrasts and reference from mapping** after any mapping edit, then select planned contrasts.
4. Choose **Full model**, **Reduced model**, or adjust the two optional random-effect checkboxes in the analysis options.
5. Choose the response transform and optional normalization reference; check the backend if needed.
6. Click **Run analysis and save XLSX**, choose a new output path, and review the reported model status.

The suggested output name is `<input-name>_zetastats-v5-analysis.xlsx` for full, `<input-name>_zetastats-v5-reduced-analysis.xlsx` for reduced, or `<input-name>_zetastats-v5-custom-analysis.xlsx` for custom settings. Background execution displays status and elapsed time.

Run from the command line:

```bash
python hierarchical_statistical_analyzer.py --input data.csv --output results.xlsx
python hierarchical_statistical_analyzer.py --input data.xlsx --sheet Measurements --transform log2 --normalize-to WT 0J
python hierarchical_statistical_analyzer.py --input data.csv --random-effects reduced
python hierarchical_statistical_analyzer.py --input data.csv --drop-replicate-treatment
python hierarchical_statistical_analyzer.py --help
```

The CLI defaults to `--random-effects full`. The reduced preset excludes both optional terms; either drop switch excludes its named term from the chosen preset. The same switches can initialize the GUI when `--input` is omitted. In the Python `analyze()` API, use `include_replicate_condition` and `include_replicate_treatment` booleans (both default to `True`).

The CLI infers mapping from column names and selects all pairwise contrasts; use the GUI or Python `analyze()` API for a restricted planned family. Exit status is 0 for an OK fit, 2 when a workbook is saved with inference withheld, and 1 for a handled input/backend error. The input file cannot be overwritten by the output.

## Results

Every exported workbook contains exactly these sheets:

| Sheet | Contents |
| --- | --- |
| `1_Descriptive_stats` | Equally weighted biological replicate means and separate within-sample summaries; raw, normalized, and transformed scales as applicable |
| `2_Variance_TypeIII` | Type III tests, random/residual variance components, conditional residual summaries, model status, selected random-effects structure and exact formula, settings, warnings, and R/package versions |
| `3_Contrasts` | Model estimated means, B-minus-A estimates, SE, individual 95% CI, Satterthwaite t/df, raw and Holm p-values, effect summaries, sample counts, status, and significance symbols |

Main descriptive SD and SEM summarize biological replicate means. Within-sample SD and SEM describe measurement variability, not biological precision. Geometric means are undefined when their input values are nonpositive.

Positive contrast estimates mean B is higher than A on the fitted response scale. For log-transformed models, estimates and CI endpoints are also back-transformed as B/A fold changes. These are model-based ratios on the log scale, not arithmetic mean ratios.

Approximate partial eta-squared is computed as `F*NumDF/(F*NumDF+DenDF)` for non-intercept Type III tests and `t^2/(t^2+DenDF)` for contrasts. It is not classical total eta-squared or total explained variance. The intercept test compares the equal-weight grand mean with zero; it is not an interaction test. A nonsignificant interaction does not establish the absence of interaction.

`Symbol` is the final contrast column and depends only on `p_holm`: `****` ≤ 0.0001, `***` ≤ 0.001, `**` ≤ 0.01, `*` ≤ 0.05, otherwise `ns`. Missing adjusted p-values give blank symbols. Review `Status` before interpreting any result; fitted estimates may still be present when inference is withheld.

## Citation

Use [`CITATION.cff`](CITATION.cff) and identify the exact software version used. The DOI badge links to the Zenodo concept record for **all versions**; use the version-specific archived DOI when available.

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE).
