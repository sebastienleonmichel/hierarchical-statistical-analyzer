
"""
Hierarchical Statistical Analyzer
=================================

GUI-based framework for:
- stratified nonparametric testing
- linear mixed-effects modeling
- replicate-aware statistical analysis
- skewed biological datasets

Supports:
- blocked experimental designs
- hierarchical inference
- multiple testing correction
- contextual data mapping

Author: Sébastien Terreau
Year: 2026
"""


import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd

from scipy.stats import rankdata

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


# ============================================================
# RESHAPE
# ============================================================

def reshape_long(df_wide, col_to_group, col_to_repl):

    rows = []

    for col in df_wide.columns:

        grp = col_to_group[col]
        rep = col_to_repl[col]

        vals = pd.to_numeric(
            df_wide[col],
            errors="coerce"
        ).dropna()

        vals = vals[vals > 0]

        for v in vals:

            rows.append({

                "value": float(v),
                "group": grp,
                "replicate_id": rep,
                "column": col

            })

    return pd.DataFrame(rows)


# ============================================================
# WILCOXON UTILITIES
# ============================================================

def mannwhitney_u(x, gA_mask):

    r = rankdata(x, method="average")

    n1 = int(gA_mask.sum())

    U = r[gA_mask].sum() - n1 * (n1 + 1) / 2

    return float(U)


def van_elteren_test(
    df_long,
    groupA,
    groupB,
    n_perm=10000,
    seed=0
):

    rng = np.random.default_rng(seed)

    observed = 0.0

    for _, sub in df_long.groupby("replicate_id"):

        sub = sub[
            sub["group"].isin([groupA, groupB])
        ]

        if len(sub) == 0:
            continue

        x = sub["value"].to_numpy()

        gA = (
            sub["group"].to_numpy() == groupA
        )

        if gA.sum() == 0:
            continue

        observed += mannwhitney_u(x, gA)

    perms = []

    for _ in range(n_perm):

        total = 0.0

        for _, sub in df_long.groupby("replicate_id"):

            sub = sub[
                sub["group"].isin([groupA, groupB])
            ]

            if len(sub) == 0:
                continue

            x = sub["value"].to_numpy()

            nA = np.sum(
                sub["group"].to_numpy() == groupA
            )

            idx = np.arange(len(x))

            rng.shuffle(idx)

            mask = np.zeros(
                len(x),
                dtype=bool
            )

            mask[idx[:nA]] = True

            total += mannwhitney_u(x, mask)

        perms.append(total)

    perms = np.array(perms)

    center = perms.mean()

    p = (
        np.sum(
            np.abs(perms - center)
            >= np.abs(observed - center)
        ) + 1
    ) / (n_perm + 1)

    return observed, p


# ============================================================
# LMM
# ============================================================

def run_lmm(
    df_long,
    contrasts,
    log_transform=False
):

    df = df_long.copy()

    if log_transform:

        df["response"] = np.log(df["value"])

    else:

        df["response"] = df["value"]

    rows = []

    raw_p = []

    for A, B in contrasts:

        sub = df[
            df["group"].isin([A, B])
        ].copy()

        # --------------------------------------------
        # REQUIRE MULTIPLE REPLICATES
        # --------------------------------------------

        if sub["replicate_id"].nunique() < 2:

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    "LMM failed: <2 replicate levels"

            })

            raw_p.append(np.nan)

            continue

        sub["group"] = (
            sub["group"]
            .astype("category")
        )

        try:

            model = smf.mixedlm(
                "response ~ group",
                data=sub,
                groups=sub["replicate_id"]
            )

            fit = model.fit(reml=False)

            # ----------------------------------------
            # ROBUST EXTRACTION
            # ----------------------------------------

            beta = float(
                fit.params.iloc[1]
            )

            pval = float(
                fit.pvalues.iloc[1]
            )

            raw_p.append(pval)

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    "Linear Mixed Model",


                "p_raw":
                    pval

            })

        except Exception as e:

            rows.append({

                "Contrast":
                    f"{A} vs {B}",

                "Method":
                    f"LMM failed: {e}"

            })

            raw_p.append(np.nan)

    # ========================================================
    # MULTIPLE TESTING CORRECTIONS
    # ========================================================

    raw_p_array = np.array(raw_p)

    valid = np.isfinite(raw_p_array)

    holm_adj = np.full(len(raw_p), np.nan)
    bonf_adj = np.full(len(raw_p), np.nan)
    fdr_adj = np.full(len(raw_p), np.nan)

    if np.sum(valid) > 0:

        _, p_holm, _, _ = multipletests(
            raw_p_array[valid],
            method="holm"
        )

        _, p_bonf, _, _ = multipletests(
            raw_p_array[valid],
            method="bonferroni"
        )

        _, p_fdr, _, _ = multipletests(
            raw_p_array[valid],
            method="fdr_bh"
        )

        holm_adj[valid] = p_holm
        bonf_adj[valid] = p_bonf
        fdr_adj[valid] = p_fdr

    for i in range(len(rows)):

        rows[i]["p_holm"] = holm_adj[i]
        rows[i]["p_bonferroni"] = bonf_adj[i]
        rows[i]["p_fdr_bh"] = fdr_adj[i]

    return pd.DataFrame(rows)


# ============================================================
# GUI
# ============================================================

class HierarchicalStatisticalAnalyzerGUI(tk.Tk):

    def __init__(self):

        super().__init__()

        self.title(
            "Hierarchical Statistical Analyzer"
        )

        self.geometry("1300x900")

        self.df = None

        self.group_vars = {}
        self.repl_vars = {}
        self.contrast_vars = {}

        self.engine_var = tk.StringVar(
            value="Both"
        )

        self.log_var = tk.BooleanVar(
            value=False
        )

        self._build_gui()

    def _build_gui(self):

        nb = ttk.Notebook(self)

        nb.pack(fill=tk.BOTH, expand=True)

        self.tab_input = ttk.Frame(nb)
        self.tab_mapping = ttk.Frame(nb)
        self.tab_contrast = ttk.Frame(nb)
        self.tab_options = ttk.Frame(nb)
        self.tab_run = ttk.Frame(nb)

        nb.add(self.tab_input, text="1) Input")
        nb.add(self.tab_mapping, text="2) Mapping")
        nb.add(self.tab_contrast, text="3) Contrasts")
        nb.add(self.tab_options, text="4) Options")
        nb.add(self.tab_run, text="5) Run")

        self.build_input()
        self.build_mapping()
        self.build_contrasts()
        self.build_options()
        self.build_run()

    def build_input(self):

        ttk.Button(
            self.tab_input,
            text="Open CSV",
            command=self.open_csv
        ).pack(padx=20, pady=20)

        self.columns_box = tk.Listbox(
            self.tab_input,
            width=100,
            height=25
        )

        self.columns_box.pack(
            padx=20,
            pady=20
        )

    def open_csv(self):

        path = filedialog.askopenfilename(
            filetypes=[("CSV", "*.csv")]
        )

        if not path:
            return

        self.df = pd.read_csv(path)

        self.columns_box.delete(0, tk.END)

        for c in self.df.columns:

            self.columns_box.insert(tk.END, c)

        self.refresh_mapping()
        self.refresh_contrasts()

    def build_mapping(self):

        self.mapping_frame = ttk.Frame(
            self.tab_mapping
        )

        self.mapping_frame.pack(
            fill=tk.BOTH,
            expand=True,
            padx=20,
            pady=20
        )

    def refresh_mapping(self):

        for w in self.mapping_frame.winfo_children():
            w.destroy()

        if self.df is None:
            return

        for i, col in enumerate(self.df.columns):

            ttk.Label(
                self.mapping_frame,
                text=col,
                width=40
            ).grid(row=i, column=0)

            g = tk.StringVar(
                value=col.split("_")[0]
            )

            r = tk.StringVar(
                value=col.split("_")[-1]
            )

            self.group_vars[col] = g
            self.repl_vars[col] = r

            ttk.Entry(
                self.mapping_frame,
                textvariable=g,
                width=20
            ).grid(row=i, column=1)

            ttk.Entry(
                self.mapping_frame,
                textvariable=r,
                width=20
            ).grid(row=i, column=2)

    def build_contrasts(self):

        self.contrast_frame = ttk.Frame(
            self.tab_contrast
        )

        self.contrast_frame.pack(
            fill=tk.BOTH,
            expand=True,
            padx=20,
            pady=20
        )

    def refresh_contrasts(self):

        for w in self.contrast_frame.winfo_children():
            w.destroy()

        if self.df is None:
            return

        groups = sorted(set([
            c.split("_")[0]
            for c in self.df.columns
        ]))

        row = 0

        for i in range(len(groups)):

            for j in range(i + 1, len(groups)):

                A = groups[i]
                B = groups[j]

                var = tk.BooleanVar(
                    value=True
                )

                self.contrast_vars[(A, B)] = var

                ttk.Checkbutton(
                    self.contrast_frame,
                    text=f"{A} vs {B}",
                    variable=var
                ).grid(row=row, column=0, sticky="w")

                row += 1

    def build_options(self):

        frm = self.tab_options

        ttk.Label(
            frm,
            text="Statistical engine"
        ).pack(anchor="w", padx=20, pady=10)

        ttk.Combobox(
            frm,
            textvariable=self.engine_var,
            values=[
                "Stratified Wilcoxon",
                "Linear Mixed Model",
                "Both"
            ],
            state="readonly"
        ).pack(anchor="w", padx=20)

        ttk.Checkbutton(
            frm,
            text="Log-transform values",
            variable=self.log_var
        ).pack(anchor="w", padx=20, pady=20)

    def build_run(self):

        ttk.Button(
            self.tab_run,
            text="Run Analysis",
            command=self.run_analysis
        ).pack(anchor="w", padx=20, pady=20)

        self.output = tk.Text(
            self.tab_run,
            width=180,
            height=40
        )

        self.output.pack(
            padx=20,
            pady=20
        )

    def run_analysis(self):

        col_to_group = {
            c: self.group_vars[c].get()
            for c in self.df.columns
        }

        col_to_repl = {
            c: self.repl_vars[c].get()
            for c in self.df.columns
        }

        df_long = reshape_long(
            self.df,
            col_to_group,
            col_to_repl
        )

        contrasts = []

        for (A, B), var in self.contrast_vars.items():

            if var.get():

                contrasts.append((A, B))

        all_results = []

        engine = self.engine_var.get()

        # ====================================================
        # STRATIFIED WILCOXON
        # ====================================================

        if engine in [
            "Stratified Wilcoxon",
            "Both"
        ]:

            rows = []

            raw_p = []

            for A, B in contrasts:

                U, p = van_elteren_test(
                    df_long,
                    A,
                    B
                )

                raw_p.append(p)

                rows.append({

                    "Contrast":
                        f"{A} vs {B}",

                    "Method":
                        "Stratified Wilcoxon",

                    "p_raw":
                        p

                })

            raw_p_array = np.array(raw_p)

            valid = np.isfinite(raw_p_array)

            holm_adj = np.full(len(raw_p), np.nan)
            bonf_adj = np.full(len(raw_p), np.nan)
            fdr_adj = np.full(len(raw_p), np.nan)

            if np.sum(valid) > 0:

                _, p_holm, _, _ = multipletests(
                    raw_p_array[valid],
                    method="holm"
                )

                _, p_bonf, _, _ = multipletests(
                    raw_p_array[valid],
                    method="bonferroni"
                )

                _, p_fdr, _, _ = multipletests(
                    raw_p_array[valid],
                    method="fdr_bh"
                )

                holm_adj[valid] = p_holm
                bonf_adj[valid] = p_bonf
                fdr_adj[valid] = p_fdr

            for i in range(len(rows)):

                rows[i]["p_holm"] = holm_adj[i]
                rows[i]["p_bonferroni"] = bonf_adj[i]
                rows[i]["p_fdr_bh"] = fdr_adj[i]

            wilcox_df = pd.DataFrame(rows)

            all_results.append(wilcox_df)

        # ====================================================
        # LMM
        # ====================================================

        if engine in [
            "Linear Mixed Model",
            "Both"
        ]:

            lmm_df = run_lmm(
                df_long,
                contrasts,
                log_transform=self.log_var.get()
            )

            all_results.append(lmm_df)

        results = pd.concat(
            all_results,
            ignore_index=True
        )

        results.to_csv(
            "hierarchical_analysis_results.csv",
            index=False
        )

        self.output.delete("1.0", tk.END)

        self.output.insert(
            tk.END,
            results.to_string(index=False)
        )

        messagebox.showinfo(
            "Done",
            "Analysis completed"
        )


if __name__ == "__main__":

    app = HierarchicalStatisticalAnalyzerGUI()

    app.mainloop()
