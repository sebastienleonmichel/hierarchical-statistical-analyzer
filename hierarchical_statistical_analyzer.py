"""ZetaStats v5.0.0 — replicate-aware factorial linear mixed models.

Author: Sébastien Terreau, 2026.
Input: wide CSV/XLSX, one sample per Condition_treatment_replicateID column.
Python provides the GUI/data handling. R's lme4/lmerTest provide REML fitting,
Type III F tests and Satterthwaite contrast tests from ONE global model.

Default two-factor model:
  response ~ condition * treatment + (1|replicate_id)
           + (1|replicate_id:condition) + (1|replicate_id:treatment)
           + (1|replicate_id:condition:treatment)

Condition and the overall treatment effect are fixed. Random replicate ×
treatment deviations let treatment effects vary between experiments, with
an exchangeable categorical covariance structure. Replicate × condition and
sample deviations represent the other pairwise/three-factor random effects.
Cells remain nested in samples; cell counts are NOT biological replicate n.
No automatic simpler-model fallback is made for invalid or singular fits.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
from pathlib import Path
import queue
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import pandas as pd

APP_VERSION = "v5.0.0"
LMM_RESPONSE_TRANSFORMS = ("None", "log2", "ln", "log10")
SHEET_NAMES = ("1_Descriptive_stats", "2_Variance_TypeIII", "3_Contrasts")
REQUIRED_R_PACKAGES = ("lme4", "lmerTest", "jsonlite")


R_ENGINE = r'''
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Expected a configuration JSON path")
suppressPackageStartupMessages(library(jsonlite))
cfg <- jsonlite::fromJSON(args[1], simplifyVector=FALSE)
out <- cfg$output
tryCatch({
  suppressPackageStartupMessages(library(lme4))
  suppressPackageStartupMessages(library(lmerTest))
  dat <- read.csv(cfg$data, stringsAsFactors=FALSE, check.names=FALSE,
                  na.strings=character(), fileEncoding="UTF-8")
  dat$condition <- factor(dat$condition, levels=unlist(cfg$conditions))
  dat$treatment <- factor(dat$treatment, levels=unlist(cfg$treatments))
  dat$replicate_id <- factor(dat$replicate_id)
  dat$sample_id <- factor(dat$sample_id)
  nC <- nlevels(dat$condition); nT <- nlevels(dat$treatment)
  if (nC > 1L) contrasts(dat$condition) <- contr.sum(nC)
  if (nT > 1L) contrasts(dat$treatment) <- contr.sum(nT)

  fixed <- if (nC > 1L && nT > 1L) "condition * treatment" else
           if (nC > 1L) "condition" else "treatment"
  random <- c("(1|replicate_id)")
  if (nC > 1L && nT > 1L) {
    random <- c(random, "(1|replicate_id:condition)",
                "(1|replicate_id:treatment)")
  }
  random <- c(random, "(1|sample_id)")
  formula_text <- paste("response ~", fixed, "+", paste(random, collapse=" + "))
  fixed_formula <- as.formula(paste("~", fixed))
  X <- model.matrix(fixed_formula, dat)
  if (qr(X)$rank != ncol(X)) stop("Fixed effects are confounded or not estimable")
  notes <- character()
  fit <- withCallingHandlers(
    lmerTest::lmer(as.formula(formula_text), data=dat, REML=TRUE,
      control=lme4::lmerControl(optimizer="bobyqa",
        optCtrl=list(maxfun=200000), check.rankX="stop.deficient")),
    warning=function(w) {notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")},
    message=function(m) {notes <<- c(notes, conditionMessage(m)); invokeRestart("muffleMessage")})
  singular <- lme4::isSingular(fit, tol=1e-4)
  conv <- unlist(fit@optinfo$conv$lme4$messages)
  opt <- unlist(fit@optinfo$conv$opt)
  reasons <- character()
  if (singular) reasons <- c(reasons, "singular / boundary random-effect variance")
  if (length(conv)) reasons <- c(reasons, conv)
  if (length(opt) && any(opt != 0)) reasons <- c(reasons, "optimizer did not converge")
  if (any(grepl("converg|Hessian|eigenvalue|unable to evaluate", notes, ignore.case=TRUE)))
    reasons <- c(reasons, notes[grepl("converg|Hessian|eigenvalue|unable to evaluate", notes, ignore.case=TRUE)])
  V <- as.matrix(vcov(fit))
  if (any(!is.finite(V)) || any(!is.finite(fixef(fit))) ||
      min(eigen(V, symmetric=TRUE, only.values=TRUE)$values) <= 0)
    reasons <- c(reasons, "invalid fixed-effect covariance")
  inference_ok <- length(reasons) == 0L
  status <- if (inference_ok) "OK" else
    paste("Inference withheld:", paste(unique(reasons), collapse="; "))
  beta <- fixef(fit)
  vc <- as.data.frame(VarCorr(fit))
  variance <- lapply(seq_len(nrow(vc)), function(i) list(
    Effect=as.character(vc$grp[i]), Variance=vc$vcov[i], SD=vc$sdcor[i],
    Variance_fraction=vc$vcov[i]/sum(vc$vcov), Status=status))
  fixed_rows <- list()
  if (inference_ok) {
    av <- anova(fit, type=3, ddf="Satterthwaite")
    L <- rep(0, length(beta)); L[which(names(beta)=="(Intercept)")] <- 1
    it <- lmerTest::contestMD(fit, L=matrix(L, nrow=1), ddf="Satterthwaite")
    av <- rbind(it[,colnames(av),drop=FALSE], av)
    rownames(av)[1] <- "Intercept"
    for (i in seq_len(nrow(av))) {
      rr <- av[i,,drop=FALSE]; effect <- rownames(av)[i]
      f <- as.numeric(rr[["F value"]]); d1 <- as.numeric(rr$NumDF); d2 <- as.numeric(rr$DenDF)
      p <- as.numeric(rr[["Pr(>F)"]])
      valid <- all(is.finite(c(f,d1,d2,p))) && d1 > 0 && d2 > 0 && f >= 0 && p >= 0 && p <= 1
      fixed_rows[[length(fixed_rows)+1L]] <- list(
        Effect=effect, Sum_sq=if(valid) as.numeric(rr[["Sum Sq"]]) else NA_real_,
        Mean_sq=if(valid) as.numeric(rr[["Mean Sq"]]) else NA_real_,
        NumDF=if(valid) d1 else NA_real_, DenDF=if(valid) d2 else NA_real_,
        F=if(valid) f else NA_real_, p_raw=if(valid) p else NA_real_,
        Eta_squared_partial=if(valid && effect != "Intercept") f*d1/(f*d1+d2) else NA_real_,
        Status=if(valid) status else "Inference withheld: invalid Type III test")
    }
  } else {
    effects <- c("Intercept", if(nC>1) "condition", if(nT>1) "treatment",
                 if(nC>1 && nT>1) "condition:treatment")
    fixed_rows <- lapply(effects, function(e) list(Effect=e, Status=status))
  }
  for (e in c("condition", "treatment", "condition:treatment")) {
    if (!(e %in% vapply(fixed_rows, function(r) r$Effect, character(1))))
      fixed_rows[[length(fixed_rows)+1L]] <- list(Effect=e,
        Status="Not applicable: this factor has only one level")
  }
  cells <- cfg$cells
  cell_design <- lapply(cells, function(cell) {
    nd <- data.frame(condition=factor(cell$condition, levels=levels(dat$condition)),
                     treatment=factor(cell$treatment, levels=levels(dat$treatment)))
    if(nC>1) contrasts(nd$condition) <- contr.sum(nC)
    if(nT>1) contrasts(nd$treatment) <- contr.sum(nT)
    row <- model.matrix(fixed_formula, nd)
    as.numeric(row[1,names(beta)])
  })
  names(cell_design) <- vapply(cells, function(cell) cell$id, character(1))
  estimated_means <- lapply(cells, function(cell) list(id=cell$id,
    EMM=as.numeric(crossprod(cell_design[[cell$id]], beta))))
  contrast_rows <- lapply(cfg$contrasts, function(cc) {
    L <- cell_design[[cc$B]] - cell_design[[cc$A]]
    rr <- list(A=cc$A, B=cc$B, Estimate_B_minus_A=as.numeric(crossprod(L,beta)),
               Status=status)
    if (inference_ok) {
      ct <- tryCatch(lmerTest::contest1D(fit, L=L, ddf="Satterthwaite",
                 confint=TRUE, level=0.95), error=function(e) e)
      if (inherits(ct,"error")) {
        rr$Status <- paste("Inference withheld:", conditionMessage(ct))
      } else {
        valid <- all(is.finite(as.numeric(ct[1,]))) && ct$df > 0 &&
                 ct[["Std. Error"]] > 0 && ct[["Pr(>|t|)"]] >= 0 && ct[["Pr(>|t|)"]] <= 1
        if (valid) {
          rr$SE <- as.numeric(ct[["Std. Error"]]); rr$DenDF <- as.numeric(ct$df)
          rr$t <- as.numeric(ct[["t value"]]); rr$p_raw <- as.numeric(ct[["Pr(>|t|)"]])
          rr$CI95_low <- as.numeric(ct[["lower"]]); rr$CI95_high <- as.numeric(ct[["upper"]])
          rr$Eta_squared_partial <- rr$t^2/(rr$t^2+rr$DenDF)
        } else rr$Status <- "Inference withheld: invalid contrast test"
      }
    }
    rr
  })
  # NA tests still count in the planned family. Never shrink multiplicity after failures.
  ps <- vapply(contrast_rows, function(r) if(is.null(r$p_raw)) NA_real_ else r$p_raw, numeric(1))
  adj <- p.adjust(ps, method="holm", n=length(ps))
  contrast_rows <- lapply(seq_along(contrast_rows), function(i) {
    rr <- contrast_rows[[i]]; rr$p_holm <- adj[i]; rr
  })
  residuals <- resid(fit)
  residual_rows <- lapply(cells, function(cell) {
    z <- residuals[dat$condition==cell$condition & dat$treatment==cell$treatment]
    list(id=cell$id, N_observations=length(z), Residual_mean=mean(z), Residual_SD=sd(z))
  })
  meta <- list(formula=formula_text, estimation="REML", df_method="Satterthwaite",
    status=status, inference_ok=inference_ok, singular=singular,
    observations=nrow(dat), replicates=nlevels(dat$replicate_id),
    samples=nlevels(dat$sample_id), warnings=as.list(unique(notes)),
    R_version=R.version.string,
    lme4_version=as.character(packageVersion("lme4")),
    lmerTest_version=as.character(packageVersion("lmerTest")),
    jsonlite_version=as.character(packageVersion("jsonlite")))
  write_json(list(meta=meta, type_iii=fixed_rows, variance=variance,
    contrasts=contrast_rows, means=estimated_means, residuals=residual_rows),
    out, auto_unbox=TRUE, pretty=TRUE, na="null", digits=17)
}, error=function(e) {
  write_json(list(error=conditionMessage(e)), out, auto_unbox=TRUE, pretty=TRUE)
  quit(status=1L)
})
'''


def parse_column_name(name):
    """Parse from the right: conditions can contain underscores."""
    parts = str(name).rsplit("_", 2)
    if len(parts) != 3 or any(not p.strip() for p in parts):
        raise ValueError(f"'{name}' must follow Condition_treatment_replicateID")
    return tuple(p.strip() for p in parts)


def read_input(path, sheet_name=0):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        # Check original headers before pandas can rename duplicate columns.
        with path.open(encoding="utf-8-sig", newline="") as fh:
            headers = next(csv.reader(fh), [])
        if len(headers) != len(set(headers)):
            raise ValueError("Duplicate column names: assign each sample a unique header.")
        df = pd.read_csv(path, keep_default_na=False, encoding="utf-8-sig")
    elif path.suffix.lower() == ".xlsx":
        from openpyxl import load_workbook
        wb = load_workbook(path, read_only=True, data_only=True)
        try:
            ws = wb.worksheets[sheet_name] if isinstance(sheet_name, int) else wb[sheet_name]
            headers = list(next(ws.iter_rows(min_row=1, max_row=1, values_only=True), []))
            used = [v for v in headers if v is not None]
            if len(used) != len(set(used)):
                raise ValueError("Duplicate column names: assign each sample a unique header.")
        finally:
            wb.close()
        df = pd.read_excel(path, sheet_name=sheet_name, engine="openpyxl", keep_default_na=False)
    else:
        raise ValueError("Select a .csv or .xlsx file.")
    if df.empty or len(df.columns) == 0:
        raise ValueError("The input contains no observations.")
    return df


def infer_mapping(df_wide):
    return {col: parse_column_name(col) for col in df_wide.columns}


def reshape_long(df_wide, mapping):
    """Mapping values are (condition, treatment, replicate); None excludes a column."""
    frames, used_samples = [], set()
    for col in df_wide.columns:
        if col not in mapping:
            raise ValueError(f"No mapping for column '{col}'.")
        if mapping[col] is None:
            continue
        if len(mapping[col]) != 3:
            raise ValueError(f"Provide condition, treatment and replicate ID for '{col}'.")
        c, t, r = (str(v).strip() for v in mapping[col])
        if not c or not t or not r:
            raise ValueError(f"Missing condition, treatment or replicate ID for '{col}'.")
        sample = (c, t, r)
        if sample in used_samples:
            raise ValueError(f"Multiple columns map to {sample}. One column must represent one sample; combine technical measurements explicitly before import.")
        used_samples.add(sample)
        series = df_wide[col]
        blank = series.isna() | series.astype(str).str.strip().eq("")
        numbers = pd.to_numeric(series, errors="coerce")
        bad = ~blank & (numbers.isna() | ~np.isfinite(numbers))
        if bad.any():
            raise ValueError(f"'{col}' has {int(bad.sum())} nonnumeric or infinite observation(s). Only empty cells are omitted.")
        values = numbers[~blank].to_numpy(dtype=float)
        if len(values) == 0:
            raise ValueError(f"'{col}' contains no numeric observations. Exclude it or provide data.")
        frames.append(pd.DataFrame(dict(value=values, condition=c, treatment=t,
            replicate_id=r, column=str(col),
            sample_id=json.dumps(sample, ensure_ascii=False))))
    if not frames:
        raise ValueError("Include at least one sample column.")
    return pd.concat(frames, ignore_index=True)


def transform_lmm_response(values, response_transform="None"):
    if response_transform not in LMM_RESPONSE_TRANSFORMS:
        raise ValueError(f"Unknown response transform: {response_transform}")
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Response values must be finite.")
    if response_transform == "None":
        return values.copy()
    if (values <= 0).any():
        raise ValueError(f"{response_transform} requires strictly positive values; found {int((values <= 0).sum())} nonpositive observation(s). No pseudocount is added.")
    return {"log2": np.log2, "ln": np.log, "log10": np.log10}[response_transform](values)


def prepare_response(df_long, response_transform="None", normalize_to=None):
    df = df_long.copy()
    divisor = 1.0
    if normalize_to is not None:
        c, t = normalize_to
        ref = df[(df.condition == c) & (df.treatment == t)]
        if ref.empty:
            raise ValueError(f"Normalization reference {(c,t)} is absent.")
        means = ref.groupby("replicate_id", sort=False).value.mean()
        divisor = float(means.mean())
        if not math.isfinite(divisor) or divisor <= 0:
            raise ValueError("The reference mean of replicate means must be finite and positive.")
    df["normalized_value"] = df.value / divisor
    df["response"] = transform_lmm_response(df.normalized_value, response_transform)
    return df, divisor


def significance_symbol(p):
    if p is None or not pd.notna(p) or not math.isfinite(float(p)):
        return ""
    for cutoff, symbol in ((0.0001,"****"),(0.001,"***"),(0.01,"**"),(0.05,"*")):
        if p <= cutoff:
            return symbol
    return "ns"


def summary_values(values):
    z = np.asarray(values, dtype=float)
    n = len(z)
    sd = float(np.std(z, ddof=1)) if n > 1 else np.nan
    geometric = float(np.exp(np.mean(np.log(z)))) if n and (z > 0).all() else np.nan
    return dict(N=n, Mean=float(np.mean(z)), SEM=sd/math.sqrt(n), SD=sd,
        Median=float(np.median(z)), Geometric_mean=geometric,
        Geometric_mean_status="OK" if math.isfinite(geometric) else "Undefined: nonpositive value(s)")


def descriptive_statistics(df, response_transform, normalized):
    scales = [("Raw input", "value")]
    if normalized:
        scales.append(("Normalized input", "normalized_value"))
    if response_transform != "None":
        scales.append((f"Model response ({response_transform})", "response"))
    rows = []
    for scale, column in scales:
        for (c,t), cell in df.groupby(["condition","treatment"], sort=True):
            replicate_means = cell.groupby("replicate_id")[column].mean()
            rows.append(dict(Scale=scale, Summary_level="Biological replicate means",
                Condition=c, Treatment=t, Replicate_ID="All", N_replicates=len(replicate_means),
                N_observations=len(cell), **summary_values(replicate_means),
                Note="Equal weight per replicate; SEM and SD describe replicate means."))
            for r, sample in cell.groupby("replicate_id", sort=True):
                rows.append(dict(Scale=scale, Summary_level="Within-sample observations",
                    Condition=c, Treatment=t, Replicate_ID=r, N_replicates=1,
                    N_observations=len(sample), **summary_values(sample[column]),
                    Note="Within-sample SEM describes measurement/cell variability, not biological precision."))
    return pd.DataFrame(rows)


def find_rscript(explicit=None):
    # A frozen macOS app must use its private engine regardless of PATH or
    # user environment overrides. An explicitly supplied CLI path is opt-in.
    if explicit is None and getattr(sys, "frozen", False):
        bundled = Path(sys.executable).resolve().parent.parent / "Resources" / "R" / "bin" / "zetastats-Rscript"
        if bundled.is_file() and os.access(bundled, os.X_OK):
            return str(bundled)
        raise RuntimeError("The statistical engine included with ZetaStats is missing. Extract the complete application again.")
    configured = explicit or os.environ.get("ZETASTATS_RSCRIPT")
    if configured and (not Path(configured).is_file() or not os.access(configured,os.X_OK)):
        raise RuntimeError(f"The configured Rscript is missing or not executable: {configured}")
    candidate = configured or shutil.which("Rscript")
    if candidate and Path(candidate).is_file() and os.access(candidate,os.X_OK):
        return str(Path(candidate))
    for path in ("/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript",
                 "/Library/Frameworks/R.framework/Resources/bin/Rscript"):
        if Path(path).is_file():
            return path
    raise RuntimeError("Rscript was not found. Install R and the R packages lme4, lmerTest and jsonlite. See README.md; alternatively set ZETASTATS_RSCRIPT.")


def check_backend(rscript=None):
    executable = find_rscript(rscript)
    result = subprocess.run([executable, "--vanilla", "-e",
        'p <- c("lme4","lmerTest","jsonlite"); missing <- p[!vapply(p, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]; if(length(missing)) {cat(paste(missing,collapse=", ")); quit(status=2)}; cat(R.version.string)'],
        capture_output=True, text=True, encoding="utf-8", timeout=45)
    if result.returncode:
        raise RuntimeError("R backend unavailable. Missing packages or R error: " + (result.stdout + result.stderr).strip())
    return executable, result.stdout.strip()


def validate_design(df):
    cells = sorted(set(zip(df.condition, df.treatment)))
    conditions, treatments = sorted(df.condition.unique()), sorted(df.treatment.unique())
    if len(cells) < 2:
        raise ValueError("Include at least two condition–treatment combinations.")
    missing = set(itertools.product(conditions,treatments)) - set(cells)
    if missing:
        raise ValueError(f"Missing condition–treatment combinations: {sorted(missing)}. The full interaction is not estimable.")
    if df.replicate_id.nunique() < 2:
        raise ValueError("At least two independent biological replicates are required.")
    per_cell = df.groupby(["condition","treatment"]).replicate_id.nunique()
    if (per_cell < 2).any():
        raise ValueError("Each condition–treatment combination needs at least two biological replicates.")
    if (df.groupby("sample_id").size() < 2).any():
        raise ValueError("Each sample needs repeated observations to separate sample variance from residual variance.")
    # Detect structurally identical/aliased random components before fitting.
    samples = df.drop_duplicates("sample_id")
    keys = [[(r,) for r in samples.replicate_id]]
    if len(conditions)>1 and len(treatments)>1:
        keys += [list(zip(samples.replicate_id,samples.condition)),
                 list(zip(samples.replicate_id,samples.treatment))]
    keys += [[(s,) for s in samples.sample_id]]
    kernels = np.column_stack([np.array([[a==b for b in key] for a in key], dtype=float).ravel() for key in keys])
    if np.linalg.matrix_rank(kernels) < len(keys):
        raise ValueError("Replicate/condition/treatment random components are confounded. Check replicate IDs and the repeated design.")
    return cells, conditions, treatments


def run_r_model(df, contrasts, rscript=None, timeout=1800, progress_callback=None):
    cells, conditions, treatments = validate_design(df)
    cell_ids = {cell: f"G{i+1:04d}" for i,cell in enumerate(cells)}
    executable, _ = check_backend(rscript)
    if progress_callback:
        progress_callback("Fitting one full LMM (REML); computing Type III and contrast tests...")
    with tempfile.TemporaryDirectory(prefix="zetastats-v5-") as temporary:
        folder = Path(temporary)
        data_path, config_path = folder/"data.csv", folder/"config.json"
        result_path, script_path = folder/"results.json", folder/"engine.R"
        df[["response","condition","treatment","replicate_id","sample_id"]].to_csv(data_path, index=False)
        config = dict(data=str(data_path), output=str(result_path), conditions=conditions,
            treatments=treatments, cells=[dict(id=cell_ids[c],condition=c[0],treatment=c[1]) for c in cells],
            contrasts=[dict(A=cell_ids[a],B=cell_ids[b]) for a,b in contrasts])
        config_path.write_text(json.dumps(config, ensure_ascii=False), encoding="utf-8")
        script_path.write_text(R_ENGINE, encoding="utf-8")
        process = subprocess.run([executable, "--vanilla", str(script_path), str(config_path)],
            capture_output=True, text=True, encoding="utf-8", timeout=timeout)
        result = json.loads(result_path.read_text(encoding="utf-8")) if result_path.exists() else {}
        if result.get("error") or process.returncode:
            raise RuntimeError(result.get("error") or (process.stderr or process.stdout)[-3000:] or "R model failed without a result.")
    return result, cell_ids


CONTRAST_COLUMNS = ["Contrast", "Condition_A", "Treatment_A", "Condition_B", "Treatment_B",
    "Method", "Response_transform", "Normalization_reference", "Normalization_divisor",
    "EMM_A", "EMM_B", "Estimate_B_minus_A", "SE", "CI95_low", "CI95_high", "t", "DenDF",
    "p_raw", "p_holm", "Eta_squared_partial", "Fold_change_B_over_A", "Fold_CI95_low",
    "Fold_CI95_high", "N_replicates_A", "N_replicates_B", "N_observations_A", "N_observations_B",
    "Holm_family_size", "Status", "Symbol"]


def cell_label(cell):
    return f"{cell[0]} | {cell[1]}"


def format_results(df, backend, cell_ids, contrasts, transform, reference, divisor):
    inverse = {v:k for k,v in cell_ids.items()}
    meta = backend["meta"]
    means = {r["id"]:r["EMM"] for r in backend["means"]}
    rows = []
    for result in backend["contrasts"]:
        a,b = inverse[result["A"]], inverse[result["B"]]
        row = {k:v for k,v in result.items() if k not in ("A","B")}
        row.update(Contrast=f"{cell_label(a)} vs {cell_label(b)}", Condition_A=a[0],Treatment_A=a[1],
            Condition_B=b[0],Treatment_B=b[1], Method="Global factorial LMM; Satterthwaite t test",
            Response_transform=transform, Normalization_reference=cell_label(reference) if reference else "None",
            Normalization_divisor=divisor, EMM_A=means[result["A"]], EMM_B=means[result["B"]],
            Holm_family_size=len(contrasts))
        for letter,cell in (("A",a),("B",b)):
            sample = df[(df.condition==cell[0]) & (df.treatment==cell[1])]
            row[f"N_replicates_{letter}"] = sample.replicate_id.nunique()
            row[f"N_observations_{letter}"] = len(sample)
        if transform != "None":
            base = {"ln":math.e,"log2":2.0,"log10":10.0}[transform]
            for source,target in (("Estimate_B_minus_A","Fold_change_B_over_A"),
                                  ("CI95_low","Fold_CI95_low"),("CI95_high","Fold_CI95_high")):
                v = row.get(source)
                if v is not None:
                    with np.errstate(over="ignore"):
                        fold = float(np.exp(float(v)*math.log(base)))
                    row[target] = fold if math.isfinite(fold) else np.nan
        row["Symbol"] = significance_symbol(row.get("p_holm"))
        rows.append(row)
    contrast_df = pd.DataFrame(rows).reindex(columns=CONTRAST_COLUMNS)
    variance_rows = []
    interpretations = {
        "Intercept":"Equal-weight grand mean across condition–treatment cells versus zero; this does not test interaction.",
        "condition":"Average condition effect across treatment levels, adjusted for the full model.",
        "treatment":"Average treatment effect across conditions, with random replicate-specific treatment deviations.",
        "condition:treatment":"Does the condition effect change with treatment? A nonsignificant result does not establish no interaction."}
    for rr in backend["type_iii"]:
        variance_rows.append(dict(Section="Type III fixed-effect tests", **rr,
            Note=interpretations.get(rr["Effect"],"")))
    random_names = {"replicate_id":"Replicate (N1/N2/N3)",
        "replicate_id:condition":"Replicate × condition", "condition:replicate_id":"Replicate × condition",
        "replicate_id:treatment":"Replicate × treatment", "treatment:replicate_id":"Replicate × treatment",
        "sample_id":"Replicate × condition × treatment (sample)" if df.condition.nunique()>1 and df.treatment.nunique()>1 else "Replicate × group (sample)",
        "Residual":"Within-sample residual"}
    for rr in backend["variance"]:
        rr = rr.copy(); rr["Effect"] = random_names.get(rr["Effect"],rr["Effect"])
        variance_rows.append(dict(Section="Random/residual variance components", **rr,
            Note="Variance estimate on the model response scale; no ordinary fixed-effect F test applies."))
    for rr in backend["residuals"]:
        rr = rr.copy(); cell = inverse[rr.pop("id")]
        variance_rows.append(dict(Section="Conditional residual diagnostics", Effect=cell_label(cell),
            **rr, Note="Descriptive residual variability; assess variance patterns after transformation."))
    notes = dict(App_version=APP_VERSION, Model_formula=meta["formula"],
        Estimation=meta["estimation"], DF_method=meta["df_method"], Model_status=meta["status"],
        Input_observations=len(df), Biological_replicates=df.replicate_id.nunique(), Biological_samples=df.sample_id.nunique(),
        Response_transform=transform, Normalization_reference=cell_label(reference) if reference else "None",
        Normalization_divisor=divisor,
        Normalization_definition="All observations divided by the arithmetic mean of the reference replicate means, BEFORE transformation; conditional on the observed divisor.",
        Random_treatment_definition="Independent exchangeable replicate × treatment deviations plus a fixed overall treatment effect; not a fully random treatment factor or an unrestricted correlated-slope covariance.",
        Interaction_definition="Condition × treatment is fixed; replicate × condition, replicate × treatment and replicate × condition × treatment are random variance components.",
        Eta_squared_definition="Approximate partial eta-squared: F*NumDF/(F*NumDF+DenDF); contrasts use t^2/(t^2+DenDF). Not classical total eta-squared or total explained variance.",
        P_value_definition="Type III and all contrast tests use the same full model and all included observations; no pairwise refits.",
        Holm_definition=f"Across all {len(contrasts)} selected contrasts; failed tests still count. Type III omnibus tests form a separate unadjusted table.",
        CI_definition="95% individual model confidence intervals; these intervals are not adjusted for multiplicity.",
        Descriptive_definition="Main rows summarize equally weighted biological replicate means. Separate within-sample rows describe observations.",
        Residual_assumption="Gaussian residuals with common variance on the fitted scale; logs can help but do not guarantee this assumption.",
        Sample_definition="One condition × treatment × replicate column is one biological sample; technical measurements may be combined only explicitly.",
        Input_columns="; ".join(df.column.unique()),
        Fit_warnings="; ".join(meta.get("warnings",[])),
        R_version=meta["R_version"], lme4_version=meta["lme4_version"],
        lmerTest_version=meta["lmerTest_version"], jsonlite_version=meta["jsonlite_version"])
    if df.replicate_id.nunique()<5:
        notes["Few_replicates"]="Fewer than five independent experiments: random-effect variance and denominator degrees of freedom may be imprecise. More cells do not add biological replicates."
    for key,value in notes.items():
        variance_rows.append(dict(Section="Model/settings",Effect=key,Note=str(value)))
    variance_df = pd.DataFrame(variance_rows).reindex(columns=["Section","Effect","Sum_sq","Mean_sq","NumDF","DenDF","F","p_raw",
        "Eta_squared_partial","Variance","SD","Variance_fraction","N_observations","Residual_mean","Residual_SD","Status","Note"])
    return variance_df, contrast_df, notes


def failed_results(df, contrasts, transform, reference, divisor, error):
    status = f"Inference withheld: {error}"
    rows = []
    for a,b in contrasts:
        rows.append(dict(Contrast=f"{cell_label(a)} vs {cell_label(b)}",Condition_A=a[0],Treatment_A=a[1],
            Condition_B=b[0],Treatment_B=b[1],Method="Global factorial LMM; unavailable",
            Response_transform=transform,Normalization_reference=cell_label(reference) if reference else "None",
            Normalization_divisor=divisor,Holm_family_size=len(contrasts),Status=status,Symbol=""))
    variance = pd.DataFrame([dict(Section="Model/settings",Effect="Model_status",Status=status,Note=str(error)),
        dict(Section="Model/settings",Effect="Response_transform",Note=transform),
        dict(Section="Model/settings",Effect="Normalization_reference",Note=cell_label(reference) if reference else "None"),
        dict(Section="Model/settings",Effect="Normalization_divisor",Note=str(divisor))])
    return variance, pd.DataFrame(rows).reindex(columns=CONTRAST_COLUMNS), dict(Model_status=status)


def analyze(df_wide, mapping=None, contrasts=None, response_transform="None", normalize_to=None,
            rscript=None, progress_callback=None, timeout=1800):
    if progress_callback:
        progress_callback("Reading sample mapping and preparing response...")
    df = reshape_long(df_wide, infer_mapping(df_wide) if mapping is None else mapping)
    cells = sorted(set(zip(df.condition,df.treatment)))
    contrasts = list(itertools.combinations(cells,2)) if contrasts is None else list(contrasts)
    if not contrasts:
        raise ValueError("Select at least one contrast.")
    seen = set()
    for a,b in contrasts:
        if a not in cells or b not in cells or a==b:
            raise ValueError(f"Invalid or outdated contrast: {a} vs {b}")
        key = frozenset((a,b))
        if key in seen:
            raise ValueError("A comparison was selected twice (including reversed direction).")
        seen.add(key)
    df, divisor = prepare_response(df, response_transform, normalize_to)
    descriptive = descriptive_statistics(df,response_transform,normalize_to is not None)
    try:
        backend,ids = run_r_model(df,contrasts,rscript,timeout,progress_callback)
        variance,results,notes = format_results(df,backend,ids,contrasts,response_transform,normalize_to,divisor)
    except (ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        variance,results,notes = failed_results(df,contrasts,response_transform,normalize_to,divisor,str(error))
    return dict(descriptive=descriptive,variance=variance,contrasts=results,metadata=notes)


def export_workbook(analysis, path):
    """Write exactly three sheets atomically; text from inputs cannot become formulas."""
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    path = Path(path)
    if path.suffix.lower() != ".xlsx":
        raise ValueError("The three-sheet output must have the .xlsx extension.")
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".xlsx",dir=path.parent,delete=False) as temp:
        temporary = Path(temp.name)
    try:
        with pd.ExcelWriter(temporary,engine="openpyxl") as writer:
            for key,sheet in zip(("descriptive","variance","contrasts"),SHEET_NAMES):
                table = analysis[key]
                table.to_excel(writer,sheet_name=sheet,index=False)
                ws = writer.sheets[sheet]
                ws.freeze_panes="A2"; ws.auto_filter.ref=ws.dimensions
                for cell in ws[1]:
                    cell.font=Font(bold=True,color="FFFFFF")
                    cell.fill=PatternFill("solid",fgColor="17365D")
                    cell.alignment=Alignment(vertical="center",wrap_text=True)
                ws.row_dimensions[1].height=32
                for column in ws.columns:
                    name = str(column[0].value)
                    ws.column_dimensions[get_column_letter(column[0].column)].width = min(45,max(14,len(name)+2))
                    for cell in column[1:]:
                        if isinstance(cell.value,str):
                            cell.data_type="s"
                        if name in ("Note","Status"):
                            cell.alignment=Alignment(wrap_text=True,vertical="top")
                            lines=max(1,math.ceil(len(str(cell.value or ""))/43))
                            ws.row_dimensions[cell.row].height=max(ws.row_dimensions[cell.row].height or 15,min(150,15*lines))
                        elif isinstance(cell.value,(float,int)) and not name.startswith("N") and name not in ("Holm_family_size",):
                            cell.number_format="0.000000"
                if "p_raw" in table:
                    for name in ("p_raw","p_holm"):
                        if name in table:
                            for cell in ws[get_column_letter(table.columns.get_loc(name)+1)][1:]:
                                cell.number_format="0.0000E+00"
        os.replace(temporary,path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


class ScrollableFrame(ttk.Frame):
    def __init__(self,parent):
        super().__init__(parent)
        self.canvas=tk.Canvas(self,highlightthickness=0)
        vertical=ttk.Scrollbar(self,orient="vertical",command=self.canvas.yview)
        horizontal=ttk.Scrollbar(self,orient="horizontal",command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=vertical.set,xscrollcommand=horizontal.set)
        self.canvas.grid(row=0,column=0,sticky="nsew")
        vertical.grid(row=0,column=1,sticky="ns"); horizontal.grid(row=1,column=0,sticky="ew")
        self.rowconfigure(0,weight=1); self.columnconfigure(0,weight=1)
        self.scrollable_frame=ttk.Frame(self.canvas)
        self.canvas.create_window((0,0),window=self.scrollable_frame,anchor="nw")
        self.scrollable_frame.bind("<Configure>",lambda event:self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<MouseWheel>",lambda event:self.canvas.yview_scroll(-1 if event.delta>0 else 1,"units"))
        self.canvas.bind("<Button-4>",lambda event:self.canvas.yview_scroll(-1,"units"))
        self.canvas.bind("<Button-5>",lambda event:self.canvas.yview_scroll(1,"units"))


class HierarchicalStatisticalAnalyzerGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"ZetaStats {APP_VERSION}"); self.geometry("1250x850"); self.minsize(920,640)
        self.df=None; self.input_path=None; self.mapping_vars={}; self.contrast_vars={}
        self.transform_var=tk.StringVar(value="None"); self.normalize_var=tk.BooleanVar(value=False)
        self.reference_var=tk.StringVar(); self.reference_cells={}
        self.status_var=tk.StringVar(value="Ready"); self.elapsed_var=tk.StringVar(value="Elapsed: 00:00")
        self.worker_queue=queue.Queue(); self.analysis_thread=None; self.started=None; self.busy=False
        self.mapping_generation=None
        self.protocol("WM_DELETE_WINDOW",self._close)
        self._build_gui()

    def _build_gui(self):
        notebook=ttk.Notebook(self); notebook.pack(fill="both",expand=True,padx=10,pady=10)
        for name,title in (("input","1) Input"),("mapping","2) Mapping"),("contrasts","3) Contrasts"),("options","4) Options"),("run","5) Run")):
            tab=ttk.Frame(notebook); notebook.add(tab,text=title); setattr(self,f"tab_{name}",tab)
        notebook.enable_traversal()
        for i in range(5):
            self.bind(f"<Command-Key-{i+1}>",lambda event,index=i:notebook.select(index))
            self.bind(f"<Control-Key-{i+1}>",lambda event,index=i:notebook.select(index))
        self.notebook=notebook
        self.build_input(); self.build_mapping(); self.build_contrasts(); self.build_options(); self.build_run()
        self.bind("<Command-o>",lambda event:self.open_data_file())
        self.bind("<Control-o>",lambda event:self.open_data_file())
        self.bind("<Command-Return>",lambda event:self.run_analysis())
        self.bind("<Control-Return>",lambda event:self.run_analysis())

    def report_callback_exception(self,exc,value,traceback):
        messagebox.showerror("Interface error",str(value))

    def build_input(self):
        self.open_button=ttk.Button(self.tab_input,text="Open CSV or XLSX",command=self.open_data_file)
        self.open_button.pack(anchor="w",padx=20,pady=20)
        ttk.Label(self.tab_input,text="One biological sample per column: Condition_treatment_replicateID\nExample: SCAIko_2J_N1. Use the same N1 label for samples from the same experiment.\nEach row is an observation within that sample. Blank cells are allowed.",justify="left").pack(anchor="w",padx=20,pady=5)
        self.file_label=ttk.Label(self.tab_input,text="No file selected",wraplength=1100)
        self.file_label.pack(anchor="w",padx=20,pady=10)
        self.columns_box=tk.Listbox(self.tab_input,height=20,width=110)
        self.columns_box.pack(fill="both",expand=True,padx=20,pady=20)

    def open_data_file(self):
        if self.busy:return
        path=filedialog.askopenfilename(title="Open data",filetypes=[("CSV / Excel","*.csv *.xlsx")])
        if not path:return
        try:
            df=read_input(path)
        except Exception as error:
            messagebox.showerror("Input error",str(error)); return
        self.df=df; self.input_path=Path(path)
        self.file_label.configure(text=f"{path}\n{len(df.columns)} columns; {len(df)} rows. Excel input uses the first worksheet.")
        self.columns_box.delete(0,tk.END)
        for col in df:self.columns_box.insert(tk.END,str(col))
        self.refresh_mapping(); self.refresh_contrasts()
        self.status_var.set("Input loaded. Check condition, treatment and replicate mapping.")
        self.notebook.select(self.tab_mapping)
        self.after_idle(self.focus_force)

    def build_mapping(self):
        ttk.Label(self.tab_mapping,text="Map condition and treatment separately. Exclude annotation columns. A replicate ID identifies an independent experiment.",wraplength=1100).pack(anchor="w",padx=20,pady=15)
        self.mapping_scroll=ScrollableFrame(self.tab_mapping); self.mapping_scroll.pack(fill="both",expand=True,padx=20,pady=10)
        self.mapping_frame=self.mapping_scroll.scrollable_frame

    def refresh_mapping(self):
        for widget in self.mapping_frame.winfo_children():widget.destroy()
        self.mapping_vars={}
        if self.df is None:return
        for j,text in enumerate(("Include","Input column","Condition","Treatment","Replicate ID")):
            ttk.Label(self.mapping_frame,text=text,font=("Helvetica",12,"bold")).grid(row=0,column=j,sticky="w",padx=8,pady=8)
        for i,col in enumerate(self.df,start=1):
            try:c,t,r=parse_column_name(col)
            except ValueError:c,t,r="","",""
            include=tk.BooleanVar(value=True)
            variables=(include,tk.StringVar(value=c),tk.StringVar(value=t),tk.StringVar(value=r))
            self.mapping_vars[col]=variables
            ttk.Checkbutton(self.mapping_frame,variable=include).grid(row=i,column=0,padx=8)
            ttk.Label(self.mapping_frame,text=str(col),width=38).grid(row=i,column=1,sticky="w",padx=8)
            for j,var in enumerate(variables[1:],start=2):
                ttk.Entry(self.mapping_frame,textvariable=var,width=24).grid(row=i,column=j,padx=8,pady=3)

    def _current_mapping(self):
        return {col:tuple(var.get().strip() for var in variables[1:]) if variables[0].get() else None
                for col,variables in self.mapping_vars.items()}

    def build_contrasts(self):
        toolbar=ttk.Frame(self.tab_contrasts); toolbar.pack(anchor="w",padx=20,pady=15)
        self.update_button=ttk.Button(toolbar,text="Update contrasts and reference from mapping",command=self.refresh_contrasts)
        self.update_button.pack(side="left",padx=(0,15))
        for text,value in (("Select all",True),("Clear all",False)):
            ttk.Button(toolbar,text=text,command=lambda v=value:self._set_contrasts(v)).pack(side="left",padx=5)
        ttk.Label(self.tab_contrasts,text="All selected B − A contrasts come from one model fitted to every included condition and treatment. Holm applies to the selected family.",wraplength=1100).pack(anchor="w",padx=20,pady=5)
        self.contrast_scroll=ScrollableFrame(self.tab_contrasts);self.contrast_scroll.pack(fill="both",expand=True,padx=20,pady=10)
        self.contrast_frame=self.contrast_scroll.scrollable_frame

    def _set_contrasts(self,value):
        if self.busy:return
        for var in self.contrast_vars.values():var.set(value)

    def refresh_contrasts(self):
        if self.busy:return
        for widget in self.contrast_frame.winfo_children():widget.destroy()
        old={key:var.get() for key,var in self.contrast_vars.items()}; self.contrast_vars={}
        mapping=self._current_mapping()
        cells=sorted({(v[0],v[1]) for v in mapping.values() if v is not None and v[0] and v[1]})
        for i,(a,b) in enumerate(itertools.combinations(cells,2)):
            var=tk.BooleanVar(value=old.get((a,b),True));self.contrast_vars[(a,b)]=var
            ttk.Checkbutton(self.contrast_frame,text=f"{cell_label(a)}  vs  {cell_label(b)}",variable=var).grid(row=i,column=0,sticky="w",padx=5,pady=3)
        # Numeric identifiers avoid ambiguous labels if user-entered names contain separators.
        self.reference_cells={f"{i+1}. {cell_label(cell)}":cell for i,cell in enumerate(cells)}
        self.reference_combo.configure(values=list(self.reference_cells))
        if self.reference_var.get() not in self.reference_cells:
            self.reference_var.set(next(iter(self.reference_cells),""))
        self.mapping_generation=mapping

    def build_options(self):
        frm=self.tab_options
        for text in (
            "Linear mixed model — one global condition × treatment model",
            "Fixed: condition, overall treatment, condition × treatment.\nRandom: replicate, replicate × condition, replicate × treatment, biological sample.",
            "Random treatment interpretation: treatment responses can vary between N1/N2/N3.\nCategorical random deviations share a variance; unrestricted slope correlations are not estimated.",
            "Type III F tests and contrast t tests use Satterthwaite degrees of freedom (R / lmerTest).",
        ):
            ttk.Label(frm,text=text,wraplength=1100,justify="left").pack(anchor="w",padx=20,pady=10)
        transforms=ttk.Frame(frm);transforms.pack(anchor="w",padx=20,pady=10)
        ttk.Label(transforms,text="Transform response:").pack(side="left",padx=(0,15))
        for name in LMM_RESPONSE_TRANSFORMS:
            ttk.Radiobutton(transforms,text=name,value=name,variable=self.transform_var).pack(side="left",padx=10)
        ttk.Label(frm,text="Logs require positive values. Normalization occurs before transformation. A transform can help variance patterns; inspect residual diagnostics.",wraplength=1100).pack(anchor="w",padx=20,pady=5)
        ttk.Checkbutton(frm,text="Normalize all observations to the mean of reference replicate means",variable=self.normalize_var).pack(anchor="w",padx=20,pady=15)
        refs=ttk.Frame(frm);refs.pack(anchor="w",padx=40,pady=5)
        ttk.Label(refs,text="Reference condition | treatment:").pack(side="left",padx=(0,10))
        self.reference_combo=ttk.Combobox(refs,textvariable=self.reference_var,state="readonly",width=55)
        self.reference_combo.pack(side="left")
        ttk.Label(frm,text="Each reference replicate has equal weight, regardless of its number of observations. The one shared divisor rescales every sample.",wraplength=1100).pack(anchor="w",padx=40,pady=5)
        ttk.Button(frm,text="Check statistical backend",command=self._check_backend).pack(anchor="w",padx=20,pady=20)
        ttk.Label(frm,text="Effect size: approximate partial eta-squared. Singular or failed fits keep descriptive results and withhold inference.",wraplength=1100).pack(anchor="w",padx=20,pady=5)

    def _check_backend(self):
        if self.busy:return
        def worker():
            try:
                exe,version=check_backend()
                location="Statistical engine included with ZetaStats." if getattr(sys,"frozen",False) else exe
                self.worker_queue.put(("backend",f"{version}\n{location}\nMixed-model tests are ready."))
            except Exception as error:self.worker_queue.put(("error",str(error)))
        self._set_busy(True);threading.Thread(target=worker,daemon=True).start();self.after(100,self._poll_worker_queue)

    def build_run(self):
        toolbar=ttk.Frame(self.tab_run);toolbar.pack(fill="x",padx=20,pady=15)
        self.run_button=ttk.Button(toolbar,text="Run analysis and save XLSX",command=self.run_analysis)
        self.run_button.pack(side="left",padx=(0,20))
        ttk.Label(toolbar,textvariable=self.elapsed_var).pack(side="right")
        ttk.Label(self.tab_run,textvariable=self.status_var,wraplength=1100).pack(anchor="w",padx=20,pady=5)
        self.progress_bar=ttk.Progressbar(self.tab_run,mode="indeterminate");self.progress_bar.pack(fill="x",padx=20,pady=10)
        output_frame=ttk.Frame(self.tab_run);output_frame.pack(fill="both",expand=True,padx=20,pady=10)
        self.output=tk.Text(output_frame,wrap="none",font=("Courier",11))
        vy=ttk.Scrollbar(output_frame,orient="vertical",command=self.output.yview)
        hx=ttk.Scrollbar(output_frame,orient="horizontal",command=self.output.xview)
        self.output.configure(yscrollcommand=vy.set,xscrollcommand=hx.set)
        self.output.grid(row=0,column=0,sticky="nsew");vy.grid(row=0,column=1,sticky="ns");hx.grid(row=1,column=0,sticky="ew")
        output_frame.rowconfigure(0,weight=1);output_frame.columnconfigure(0,weight=1)

    def _set_busy(self,busy):
        self.busy=busy
        for button in (self.open_button,self.update_button,self.run_button):button.configure(state="disabled" if busy else "normal")
        if busy:self.progress_bar.start(15)
        else:self.progress_bar.stop();self.started=None

    def run_analysis(self):
        if self.busy:return
        if self.df is None:
            messagebox.showerror("No input","Open a CSV or XLSX first.");return
        mapping=self._current_mapping()
        if mapping!=self.mapping_generation:
            messagebox.showerror("Mapping changed","Click 'Update contrasts and reference from mapping' before running.");return
        contrasts=[key for key,var in self.contrast_vars.items() if var.get()]
        reference=self.reference_cells.get(self.reference_var.get()) if self.normalize_var.get() else None
        try:
            reshape_long(self.df,mapping)
            if not contrasts:raise ValueError("Select at least one contrast.")
            if self.normalize_var.get() and reference is None:raise ValueError("Choose a normalization reference.")
        except Exception as error:
            messagebox.showerror("Analysis setup",str(error));return
        suggested=f"{self.input_path.stem}_zetastats-v5-analysis.xlsx"
        path=filedialog.asksaveasfilename(title="Save three-sheet analysis",initialdir=self.input_path.parent,
            initialfile=suggested,defaultextension=".xlsx",filetypes=[("Excel workbook","*.xlsx")])
        if not path:return
        if Path(path).resolve()==self.input_path.resolve():
            messagebox.showerror("Output path","Choose a new output file; the input cannot be overwritten.");return
        df_copy=self.df.copy();transform=self.transform_var.get()
        self.notebook.select(self.tab_run)
        self.output.delete("1.0",tk.END);self.status_var.set("Preparing analysis...")
        self.started=time.monotonic();self._set_busy(True)
        self.analysis_thread=threading.Thread(target=self._analysis_worker,args=(df_copy,mapping,contrasts,transform,reference,Path(path)),daemon=True)
        self.analysis_thread.start();self.after(100,self._poll_worker_queue)

    def _analysis_worker(self,df_wide,mapping,contrasts,transform,reference,path):
        try:
            results=analyze(df_wide,mapping,contrasts,transform,reference,
                progress_callback=lambda msg:self.worker_queue.put(("status",msg)))
            self.worker_queue.put(("status","Saving descriptive statistics, variance/Type III tests and contrasts..."))
            export_workbook(results,path);self.worker_queue.put(("done",results,path))
        except Exception as error:self.worker_queue.put(("error",str(error)))

    def _poll_worker_queue(self):
        if self.started is not None:
            elapsed=int(time.monotonic()-self.started);self.elapsed_var.set(f"Elapsed: {elapsed//60:02d}:{elapsed%60:02d}")
        try:
            while True:
                message=self.worker_queue.get_nowait();kind=message[0]
                if kind=="status":self.status_var.set(message[1])
                elif kind=="backend":
                    self._set_busy(False);messagebox.showinfo("Statistical backend",message[1]);self.after_idle(self.focus_force);return
                elif kind=="error":
                    self._set_busy(False);self.status_var.set("Analysis failed.");messagebox.showerror("Analysis error",message[1]);self.after_idle(self.focus_force);return
                elif kind=="done":
                    _,analysis,path=message;self._set_busy(False)
                    status=analysis["metadata"]["Model_status"]
                    self.status_var.set(f"Saved {path} — {status}")
                    columns=["Contrast","Estimate_B_minus_A","DenDF","p_raw","p_holm","Eta_squared_partial","Status","Symbol"]
                    self.output.insert(tk.END,f"{APP_VERSION}\n{path}\nModel status: {status}\n\n"+analysis["contrasts"][columns].to_string(index=False))
                    if status=="OK":messagebox.showinfo("Analysis saved",f"Three-sheet workbook saved:\n{path}")
                    else:messagebox.showwarning("Results saved; inference withheld",f"{status}\n\nDescriptive statistics and fit details are saved in:\n{path}")
                    self.after_idle(self.focus_force)
                    return
        except queue.Empty:pass
        if self.busy:self.after(100,self._poll_worker_queue)

    def _close(self):
        if self.busy:
            messagebox.showinfo("Analysis running","Please wait for the running analysis to finish before closing.");return
        self.destroy()


def main(argv=None):
    parser=argparse.ArgumentParser(description="ZetaStats v5.0.0: global factorial mixed models")
    parser.add_argument("--input",type=Path,help="Wide CSV/XLSX; omit to open the GUI")
    parser.add_argument("--output",type=Path,help="Three-sheet XLSX output")
    parser.add_argument("--transform",choices=LMM_RESPONSE_TRANSFORMS,default="None")
    parser.add_argument("--normalize-to",nargs=2,metavar=("CONDITION","TREATMENT"))
    parser.add_argument("--rscript",help="Path to Rscript")
    parser.add_argument("--sheet",default="0",help="Excel sheet name (default: first worksheet)")
    parser.add_argument("--check-backend",action="store_true")
    parser.add_argument("--version",action="version",version=APP_VERSION)
    args=parser.parse_args(argv)
    if args.check_backend:
        exe,version=check_backend(args.rscript);print(f"{version}\n{exe}\nRequired R packages available.");return 0
    if args.input is None:
        if args.rscript:os.environ["ZETASTATS_RSCRIPT"]=args.rscript
        app=HierarchicalStatisticalAnalyzerGUI();app.mainloop();return 0
    output=args.output or args.input.with_name(f"{args.input.stem}_zetastats-v5-analysis.xlsx")
    if output.resolve()==args.input.resolve():parser.error("Output cannot overwrite the input file.")
    df=read_input(args.input,0 if args.sheet=="0" else args.sheet)
    results=analyze(df,response_transform=args.transform,
        normalize_to=tuple(args.normalize_to) if args.normalize_to else None,rscript=args.rscript,
        progress_callback=lambda msg:print(msg,flush=True))
    export_workbook(results,output)
    status=results["metadata"]["Model_status"];print(f"Saved: {output}\n{status}")
    return 0 if status=="OK" else 2


if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        print(f"ZetaStats: {error}",file=sys.stderr);raise SystemExit(1)
