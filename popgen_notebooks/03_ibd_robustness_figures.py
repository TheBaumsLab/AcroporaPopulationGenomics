#!/usr/bin/env python3
"""
Supplementary figures for the IBD robustness analyses (response to Reviewer 2,
L. 171-173 and L. 189-191).

Figure S_pooled : pooled all-sites IBD, coloured by pair stratum, showing that the
                  pooled Mantel correlation spans a ~900 km gap in the distance
                  distribution and is carried by the between-subregion contrast.
Figure S_sens   : balanced-resampling sensitivity test (5 genets per site x 100
                  iterations) -- full-data vs. resampled FST, and the distribution of
                  Mantel r across iterations.

Run 03_ibd_pooled_all_sites.py and 03_ibd_subsampling_sensitivity.py first; this
script reads their outputs from $APAL_OUT_DIR (default <project>/output/ibd_robustness)
and writes the figures to $APAL_OUT_DIR/figures.
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# emit TrueType rather than Type 3 so PDF text stays editable in vector editors
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

PROJECT_DIR = Path(os.environ.get("APAL_PROJECT_DIR", "."))
RESULTS_DIR = Path(os.environ.get("APAL_OUT_DIR", PROJECT_DIR / "output" / "ibd_robustness"))
IBD_DIR = RESULTS_DIR
PLOTS_DIR = RESULTS_DIR / "figures"

STRATUM_STYLE = {
    "within_GreaterAntilles": ("#3B7DD8", "within Greater Antilles"),
    "within_Florida": ("#E08214", "within Florida Reef Tract"),
    "between_subregion": ("#9A9A9A", "between subregions"),
}


def save(fig, name):
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(PLOTS_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", PLOTS_DIR / f"{name}.pdf")


def fig_pooled():
    pairs = pd.read_csv(IBD_DIR / "pooled" / "ibd_pooled_pairs.tsv", sep="\t")
    strata = pd.read_csv(IBD_DIR / "pooled" / "ibd_pooled_strata.tsv", sep="\t").set_index("stratum")
    summ = pd.read_csv(IBD_DIR / "pooled" / "ibd_pooled_summary.tsv", sep="\t").iloc[0]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [2.2, 1]}
    )

    for stratum, (color, label) in STRATUM_STYLE.items():
        sub = pairs.loc[pairs["stratum"] == stratum]
        ax.scatter(sub["geo_dist_km"], sub["linearized_fst"], s=22, alpha=0.7,
                   color=color, edgecolor="none", label=f"{label} (n={len(sub)})")

    x = np.linspace(pairs["geo_dist_km"].min(), pairs["geo_dist_km"].max(), 200)
    row = strata.loc["pooled_all_sites"]
    inter = (pairs["linearized_fst"] - row["ols_slope"] * pairs["geo_dist_km"]).mean()
    ax.plot(x, row["ols_slope"] * x + inter, color="#4D4D4D", lw=1.4,
            label=f"pooled OLS ($R^2$={row['ols_r2']:.2f})")

    gap_lo = pairs.loc[pairs["stratum"] != "between_subregion", "geo_dist_km"].max()
    gap_hi = pairs.loc[pairs["stratum"] == "between_subregion", "geo_dist_km"].min()
    ax.axvspan(gap_lo, gap_hi, color="#D62728", alpha=0.07, zorder=0)
    ax.text((gap_lo + gap_hi) / 2, ax.get_ylim()[1] * 0.96,
            f"no site pairs\n{gap_lo:.0f}–{gap_hi:.0f} km",
            ha="center", va="top", fontsize=8, color="#B03030")

    ax.set_xlabel("Great-circle distance (km)")
    ax.set_ylabel(r"Linearized $F_{ST}$")
    ax.set_title("Pooled IBD across all sites")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)

    ax.text(0.98, 0.03,
            f"pooled Mantel r = {summ['mantel_r_pooled']:.2f} "
            f"(p = {summ['mantel_p_pooled']:.4f})\n"
            f"subregion split alone: r = {summ['mantel_r_subregion_split_only']:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))

    # per-stratum slopes
    order = ["within_Florida", "within_GreaterAntilles", "between_subregion", "pooled_all_sites"]
    labels = ["within FRT", "within GA", "between", "pooled"]
    colors = [STRATUM_STYLE.get(k, ("#4D4D4D",))[0] for k in order]
    vals = [strata.loc[k, "ols_slope"] * 1000 for k in order]
    ax2.bar(labels, vals, color=colors)
    for i, (v, k) in enumerate(zip(vals, order)):
        ax2.text(i, v, f" $R^2$={strata.loc[k, 'ols_r2']:.2f}", ha="center",
                 va="bottom", fontsize=8)
    ax2.set_ylabel(r"OLS slope (linearized $F_{ST}$ per 1000 km)")
    ax2.set_title("IBD slope by stratum")
    ax2.tick_params(axis="x", labelrotation=20)
    ax2.grid(alpha=0.25, axis="y")
    ax2.set_ylim(0, max(vals) * 1.25)

    fig.tight_layout()
    save(fig, "figS_ibd_pooled_all_sites")


def fig_sensitivity():
    """
    One row per habitat (Greater Antilles top, Florida Reef Tract bottom), three columns:
      1. pairwise FST, all genets vs. balanced replicates
      2. Mantel r across replicates, each panel on its own scale
      3. the same distributions on one common scale, so the difference in where they
         fall -- and how tightly -- is visible without reading the axes
    Habitat colours match the pooled figure so the two read as a set.
    """
    summ = pd.read_csv(IBD_DIR / "sensitivity" / "ibd_sensitivity_summary.tsv",
                       sep="\t").set_index("subregion")

    rows = [("greaterantilles", "GreaterAntilles", "Greater Antilles", "#3B7DD8", "#BBD0EC"),
            ("florida", "Florida", "Florida Reef Tract", "#E08214", "#F3D6B4")]

    data = {}
    for stem, key, _, _, _ in rows:
        data[key] = (
            pd.read_csv(IBD_DIR / "sensitivity" / f"ibd_sensitivity_pairs_{stem}.tsv", sep="\t"),
            pd.read_csv(IBD_DIR / "sensitivity" / f"ibd_sensitivity_iterations_{stem}.tsv", sep="\t"),
        )

    # common bins for the shared-scale column
    all_r = np.concatenate([d[1]["mantel_r"].to_numpy() for d in data.values()])
    lo = np.floor((all_r.min() - 0.05) * 20) / 20
    hi = np.ceil((all_r.max() + 0.05) * 20) / 20
    shared_bins = np.arange(lo, hi + 0.025, 0.025)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.2),
                             gridspec_kw={"width_ratios": [1.15, 1.0, 1.25]})
    shared_counts = []

    for r, (_, key, nice, color, light) in enumerate(rows):
        pairs, iters = data[key]
        s = summ.loc[key]

        # --- column 1: FST, all genets vs balanced
        ax = axes[r, 0]
        lo_f = min(pairs["fst_full_recomputed"].min(), pairs["fst_resampled_mean"].min())
        hi_f = max(pairs["fst_full_recomputed"].max(), pairs["fst_resampled_mean"].max())
        pad = 0.05 * (hi_f - lo_f)
        ax.plot([lo_f - pad, hi_f + pad], [lo_f - pad, hi_f + pad],
                color="#BBBBBB", lw=1, zorder=0)
        ax.errorbar(pairs["fst_full_recomputed"], pairs["fst_resampled_mean"],
                    yerr=[pairs["fst_resampled_mean"] - pairs["fst_resampled_q025"],
                          pairs["fst_resampled_q975"] - pairs["fst_resampled_mean"]],
                    fmt="o", ms=3.5, lw=0.6, alpha=0.7, color=color, ecolor=light)
        ax.set_xlabel(r"$F_{ST}$, all genets per site")
        ax.set_ylabel(r"$F_{ST}$, 5 genets per site" "\n" "(mean of 100 iterations)")
        ax.text(0.03, 0.97,
                f"r = {s['fst_full_vs_resampled_r']:.3f}\n"
                f"mean |Δ| = {s['fst_mean_abs_diff']:.4f}",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
        ax.grid(alpha=0.25)
        ax.locator_params(axis="both", nbins=5)

        # row label on the left-hand panel
        ax.annotate(f"{nice}\n{int(s['n_genets_total'])} genets, "
                    f"{int(s['n_genets_per_site_min'])}–{int(s['n_genets_per_site_max'])} per site",
                    xy=(-0.34, 0.5), xycoords="axes fraction", rotation=90,
                    ha="center", va="center", fontsize=11, color=color, weight="bold",
                    linespacing=1.5)

        # --- column 2: Mantel r, own scale
        ax = axes[r, 1]
        ax.hist(iters["mantel_r"], bins=20, color=color, alpha=0.8, edgecolor="white")
        ax.axvline(s["mantel_r_full"], color="#D62728", lw=1.6,
                   label=f"all genets, r = {s['mantel_r_full']:.2f}")
        ax.set_xlabel("Mantel r (balanced replicates)")
        ax.set_ylabel("iterations")
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        ax.text(0.97, 0.97,
                f"mean r = {s['mantel_r_iter_mean']:.2f}\n"
                f"95% range {s['mantel_r_iter_q025']:.2f} to {s['mantel_r_iter_q975']:.2f}\n"
                f"{s['frac_iter_mantel_p_lt_0.05'] * 100:.0f}% of replicates p < 0.05",
                transform=ax.transAxes, ha="right", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
        ax.grid(alpha=0.25, axis="y")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.38)

        # --- column 3: Mantel r, common scale (x and y shared across rows)
        ax = axes[r, 2]
        counts, _, _ = ax.hist(iters["mantel_r"], bins=shared_bins,
                               color=color, alpha=0.85, edgecolor="white", linewidth=0.4)
        shared_counts.append(counts.max())
        ax.axvline(s["mantel_r_full"], color="#D62728", lw=1.6,
                   label="all genets" if r == 0 else None)
        if r == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper left")
        ax.set_xlim(lo, hi)
        ax.set_xlabel("Mantel r (common scale)")
        ax.set_ylabel("iterations")
        ax.grid(alpha=0.25, axis="y")

    # equalise the shared-scale column so concentration vs. spread is directly comparable
    ymax = max(shared_counts) * 1.18
    for r in range(2):
        axes[r, 2].set_ylim(0, ymax)

    for ax, title in zip(axes[0], [r"Pairwise $F_{ST}$",
                                   "Mantel r — per-panel scale",
                                   "Mantel r — common scale"]):
        ax.set_title(title, fontsize=11, pad=10)

    fig.tight_layout(w_pad=2.0, h_pad=2.4)
    save(fig, "figS_ibd_balanced_resampling")


if __name__ == "__main__":
    fig_pooled()
    fig_sensitivity()
