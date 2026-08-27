#!/usr/bin/env python3
"""
Sensitivity of the IBD result to unbalanced sampling effort per site.

Reviewer 2 (L. 189-191) notes that genet counts per site are unbalanced and asks
whether pairwise FST -- and therefore the IBD signal -- is driven by that imbalance.
This script implements the test they suggest: randomly draw N genets per site,
repeat R times, and recompute pairwise FST and the Mantel/OLS statistics on each
balanced replicate.

Design decisions
----------------
* Resampling unit is the *genet* (one row per genet x site), matching the published
  analysis, which already collapses ramets.
* The locus panel is fixed to the panel used in the published analysis (per-subregion
  MAF >= 0.05 computed on the full subregion sample set; no further MAC or call-rate
  term -- see locus_mask_for). Recomputing MAF on 5 genets/site would confound "does
  sample size bias FST" with "does sample size change which loci pass the filter";
  fixing the panel isolates the effect the reviewer is asking about.
* Site coordinates, site merging and along-tract ordering are taken from the published
  run and held fixed -- only the genets contributing to each site vary.
* FST is Weir & Cockerham (1984) via scikit-allel, ratio-of-sums over loci, identical
  to the published implementation.

Inputs (paths follow 03_ibd_analysis.ipynb; set the env vars below):
  APAL_PROJECT_DIR  folder with the deposited IBD tables:
                      processed_data_03_apal_ibd_samples.tsv
                      processed_data_03a_apal_ibd_fst_greater_antilles.tsv
                      processed_data_03b_apal_ibd_fst_florida.tsv
  APAL_DATA_DIR     folder with the genotype cache mac_ld.gt_matrix.int8.npz
                    (generate it from dataset 3 with scripts/convert_vcf_to_npz.py)
Outputs (written to $APAL_OUT_DIR/sensitivity, default <project>/output/ibd_robustness):
  ibd_sensitivity_pairs_<subregion>.tsv    per-pair full vs. resampled FST
  ibd_sensitivity_iterations_<subregion>.tsv  per-iteration Mantel / OLS statistics
  ibd_sensitivity_summary.tsv              one row per subregion, for the response letter
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import allel
import numpy as np
import pandas as pd
from scipy import stats

PROJECT_DIR = Path(os.environ.get("APAL_PROJECT_DIR", "."))
DATA_DIR = Path(os.environ.get("APAL_DATA_DIR", PROJECT_DIR / "datasets"))
RESULTS_DIR = Path(os.environ.get("APAL_OUT_DIR", PROJECT_DIR / "output" / "ibd_robustness"))

NPZ_PATH = DATA_DIR / "mac_ld.gt_matrix.int8.npz"
SAMPLES_TSV = PROJECT_DIR / "processed_data_03_apal_ibd_samples.tsv"
FST_GA_TSV = PROJECT_DIR / "processed_data_03a_apal_ibd_fst_greater_antilles.tsv"
FST_FL_TSV = PROJECT_DIR / "processed_data_03b_apal_ibd_fst_florida.tsv"
OUT_DIR = RESULTS_DIR / "sensitivity"

# published pairwise tables, used both for validation and for fixed geography
PUBLISHED = {
    "GreaterAntilles": (FST_GA_TSV, "geo_dist_km"),
    # NB: the published Florida statistics use great-circle distance, not along-tract
    # distance; --florida-along-tract switches to the latter for comparison.
    "Florida": (FST_FL_TSV, "geo_dist_km"),
}

MAF_THRESH = 0.05  # regional MAF, matching REGIONAL_MAF_MIN in the published notebook


# --------------------------------------------------------------------------- FST


def dosage_to_genotype_array(gn):
    """Dosage matrix (variants x samples, values -1/0/1/2) -> allel.GenotypeArray."""
    n_variants, n_samples = gn.shape
    gt = np.full((n_variants, n_samples, 2), -1, dtype=np.int8)
    gt[gn == 0] = [0, 0]
    gt[gn == 1] = [0, 1]
    gt[gn == 2] = [1, 1]
    return allel.GenotypeArray(gt)


def fst_pair_with_nloci(gta, idx1, idx2):
    with np.errstate(divide="ignore", invalid="ignore"):
        a, b, c = allel.weir_cockerham_fst(gta, [idx1, idx2])
    num = a.sum(axis=1)
    den = (a + b + c).sum(axis=1)
    valid = np.isfinite(num) & np.isfinite(den) & (den > 0)
    if valid.sum() == 0:
        return np.nan, 0
    return num[valid].sum() / den[valid].sum(), int(valid.sum())


def locus_mask_for(G_sub):
    """Per-subregion locus filter, exactly as in the published analysis.

    Cells 34 and 61 of popgen_notebooks/03_ibd_analysis.ipynb apply a regional MAF
    filter and nothing else:

        _maf_mask = (_maf >= REGIONAL_MAF_MIN) & (_total > 0)

    MAF is computed by allel.count_alleles(), i.e. over *called* alleles only, so a
    locus with no calls in the subregion is excluded by the (_total > 0) term.
    Reproduces the deposited ibd_fst_*.tsv to 1e-16 (3,185/5,555 loci for
    Florida, 3,162/5,555 for the Greater Antilles).

    NB: an earlier version of this function also required MAC > 3 and call rate >= 0.9,
    copied from the exploratory notebook 03_apal_ibd_simple.ipynb. Those terms are not
    in the published path; they dropped a further 15 loci in Florida and 58 in the
    Greater Antilles and shifted pairwise FST by a median of 5e-4 / 1.4e-3.
    """
    called = G_sub >= 0
    n_called = called.sum(axis=1)
    n_chrom = 2 * n_called
    alt_ac = np.where(called, G_sub, 0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        alt_af = np.divide(
            alt_ac, n_chrom, out=np.full(alt_ac.shape, np.nan, dtype=float), where=n_chrom > 0
        )
    maf = np.minimum(alt_af, 1 - alt_af)
    return (maf >= MAF_THRESH) & (n_chrom > 0)


# ------------------------------------------------------------------------ Mantel


def mantel_test_pearson(x_mat, y_mat, n_perm=10_000, rng=None):
    """Pearson Mantel test on the upper triangle; matches the published helper."""
    rng = np.random.default_rng(42) if rng is None else rng
    x = np.asarray(x_mat, dtype=float)
    y = np.asarray(y_mat, dtype=float)
    iu = np.triu_indices_from(x, k=1)
    x_vec, y_vec = x[iu], y[iu]
    valid = np.isfinite(x_vec) & np.isfinite(y_vec)
    x_vec = x_vec[valid]
    if valid.sum() < 2:
        return np.nan, np.nan
    r_obs = stats.pearsonr(x_vec, y_vec[valid]).statistic
    n = x.shape[0]
    perm_rs = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(n)
        perm_rs[i] = stats.pearsonr(x_vec, y[perm][:, perm][iu][valid]).statistic
    p_perm = (np.sum(np.abs(perm_rs) >= np.abs(r_obs)) + 1) / (n_perm + 1)
    return r_obs, p_perm


def square_from_long(df, value_col, sites):
    mat = pd.DataFrame(np.nan, index=sites, columns=sites, dtype=float)
    for _, row in df.iterrows():
        mat.loc[row["site1"], row["site2"]] = row[value_col]
        mat.loc[row["site2"], row["site1"]] = row[value_col]
    np.fill_diagonal(mat.values, 0.0)
    return mat


def ibd_stats(pair_df, dist_col, sites, n_perm, rng):
    x = pair_df[dist_col].to_numpy(dtype=float)
    y = pair_df["linearized_fst"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    ols = stats.linregress(x[valid], y[valid])
    r, p = mantel_test_pearson(
        square_from_long(pair_df, dist_col, sites),
        square_from_long(pair_df, "linearized_fst", sites),
        n_perm=n_perm,
        rng=rng,
    )
    return {
        "n_pairs": int(valid.sum()),
        "slope": ols.slope,
        "intercept": ols.intercept,
        "ols_r2": ols.rvalue ** 2,
        "ols_p": ols.pvalue,
        "mantel_r": r,
        "mantel_p": p,
    }


# -------------------------------------------------------------------------- main


def run_subregion(subregion, samples, G, sample_ids, n_per_site, n_iter, n_perm,
                  n_perm_final, seed):
    pub_path, dist_col = PUBLISHED[subregion]
    pub = pd.read_csv(pub_path, sep="\t")

    sub = samples.loc[samples["subregion_ibd"] == subregion].copy()
    # one row per genet x site, as in the published analysis
    sub = sub.drop_duplicates(["mlg", "site_merged"]).reset_index(drop=True)

    sites = sorted(set(pub["site1"]).union(pub["site2"]))
    sub = sub.loc[sub["site_merged"].isin(sites)].reset_index(drop=True)

    id_to_col = pd.Series(np.arange(len(sample_ids)), index=sample_ids)
    sub["col"] = sub["affymetrix_id"].map(id_to_col)
    missing = sub["col"].isna().sum()
    if missing:
        raise SystemExit(f"{subregion}: {missing} samples not present in {NPZ_PATH.name}")
    sub["col"] = sub["col"].astype(int)

    cols_all = sub["col"].to_numpy()
    mask = locus_mask_for(G[:, cols_all])
    G_filt = G[mask][:, cols_all]
    print(f"[{subregion}] {len(sub)} genets, {len(sites)} sites, "
          f"{int(mask.sum())}/{G.shape[0]} loci pass filters")

    site_rows = {s: np.where(sub["site_merged"].to_numpy() == s)[0] for s in sites}
    counts = {s: len(v) for s, v in site_rows.items()}
    print(f"[{subregion}] genets per site: min {min(counts.values())}, "
          f"max {max(counts.values())}, median {int(np.median(list(counts.values())))}")

    # ---- full-data FST (validation against the published table)
    gta_full = dosage_to_genotype_array(G_filt.astype(np.int16))
    full_rows = []
    for _, row in pub.iterrows():
        f, nl = fst_pair_with_nloci(gta_full, site_rows[row["site1"]], site_rows[row["site2"]])
        full_rows.append({"site1": row["site1"], "site2": row["site2"],
                          "fst_full": f, "n_loci_full": nl})
    full = pd.DataFrame(full_rows)
    check = pub.merge(full, on=["site1", "site2"])
    delta = (check["fst"] - check["fst_full"]).abs()
    print(f"[{subregion}] reproduction of published FST: max |delta| = {delta.max():.3g}, "
          f"r = {stats.pearsonr(check['fst'], check['fst_full']).statistic:.6f}")

    # ---- balanced resampling
    rng = np.random.default_rng(seed)
    pair_keys = list(zip(pub["site1"], pub["site2"]))
    fst_iters = np.full((n_iter, len(pair_keys)), np.nan)
    iter_stats = []

    for it in range(n_iter):
        picked = {s: rng.choice(site_rows[s], size=n_per_site, replace=False) for s in sites}
        keep = np.concatenate([picked[s] for s in sites])
        keep.sort()
        remap = {old: new for new, old in enumerate(keep)}
        gta_it = dosage_to_genotype_array(G_filt[:, keep].astype(np.int16))

        vals = []
        for s1, s2 in pair_keys:
            f, _ = fst_pair_with_nloci(
                gta_it,
                np.array([remap[i] for i in picked[s1]]),
                np.array([remap[i] for i in picked[s2]]),
            )
            vals.append(f)
        fst_iters[it] = vals

        it_df = pd.DataFrame(
            {"site1": [k[0] for k in pair_keys], "site2": [k[1] for k in pair_keys],
             "fst": vals, dist_col: pub[dist_col].to_numpy()}
        )
        it_df["linearized_fst"] = it_df["fst"] / (1 - it_df["fst"])
        st = ibd_stats(it_df, dist_col, sites, n_perm, rng)
        st["iteration"] = it
        iter_stats.append(st)
        if (it + 1) % 10 == 0:
            print(f"[{subregion}]   iteration {it + 1}/{n_iter}")

    mean_fst = np.nanmean(fst_iters, axis=0)
    pairs_out = pub[["site1", "site2", "fst", dist_col]].copy()
    pairs_out = pairs_out.rename(columns={"fst": "fst_published"})
    pairs_out["fst_full_recomputed"] = full["fst_full"].to_numpy()
    pairs_out["fst_resampled_mean"] = mean_fst
    pairs_out["fst_resampled_sd"] = np.nanstd(fst_iters, axis=0, ddof=1)
    pairs_out["fst_resampled_q025"] = np.nanpercentile(fst_iters, 2.5, axis=0)
    pairs_out["fst_resampled_q975"] = np.nanpercentile(fst_iters, 97.5, axis=0)
    pairs_out["n_genets_site1"] = [counts[s] for s in pairs_out["site1"]]
    pairs_out["n_genets_site2"] = [counts[s] for s in pairs_out["site2"]]

    # ---- IBD on the averaged FST matrix (what the reviewer literally proposes)
    avg_df = pairs_out[["site1", "site2", dist_col]].copy()
    avg_df["fst"] = mean_fst
    avg_df["linearized_fst"] = avg_df["fst"] / (1 - avg_df["fst"])
    avg_stats = ibd_stats(avg_df, dist_col, sites, n_perm_final, rng)

    full_df = pairs_out[["site1", "site2", dist_col]].copy()
    full_df["fst"] = pairs_out["fst_full_recomputed"]
    full_df["linearized_fst"] = full_df["fst"] / (1 - full_df["fst"])
    full_stats = ibd_stats(full_df, dist_col, sites, n_perm_final, rng)

    iters = pd.DataFrame(iter_stats)
    ok = np.isfinite(pairs_out["fst_full_recomputed"]) & np.isfinite(mean_fst)
    summary = {
        "subregion": subregion,
        "distance_metric": dist_col,
        "n_sites": len(sites),
        "n_pairs": len(pair_keys),
        "n_genets_total": len(sub),
        "n_genets_per_site_min": min(counts.values()),
        "n_genets_per_site_max": max(counts.values()),
        "n_loci": int(mask.sum()),
        "n_per_site_resampled": n_per_site,
        "n_iterations": n_iter,
        "fst_full_mean": float(np.nanmean(pairs_out["fst_full_recomputed"])),
        "fst_resampled_mean": float(np.nanmean(mean_fst)),
        "fst_full_vs_resampled_r": float(
            stats.pearsonr(pairs_out["fst_full_recomputed"][ok], mean_fst[ok]).statistic
        ),
        "fst_mean_abs_diff": float(np.nanmean(np.abs(pairs_out["fst_full_recomputed"] - mean_fst))),
        "mantel_r_full": full_stats["mantel_r"],
        "mantel_p_full": full_stats["mantel_p"],
        "ols_r2_full": full_stats["ols_r2"],
        "mantel_r_avgfst": avg_stats["mantel_r"],
        "mantel_p_avgfst": avg_stats["mantel_p"],
        "ols_r2_avgfst": avg_stats["ols_r2"],
        "mantel_r_iter_mean": float(iters["mantel_r"].mean()),
        "mantel_r_iter_sd": float(iters["mantel_r"].std(ddof=1)),
        "mantel_r_iter_q025": float(iters["mantel_r"].quantile(0.025)),
        "mantel_r_iter_q975": float(iters["mantel_r"].quantile(0.975)),
        "frac_iter_mantel_p_lt_0.05": float((iters["mantel_p"] < 0.05).mean()),
        "ols_r2_iter_mean": float(iters["ols_r2"].mean()),
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tag = subregion.lower()
    pairs_out.to_csv(OUT_DIR / f"ibd_sensitivity_pairs_{tag}.tsv", sep="\t", index=False)
    iters.to_csv(OUT_DIR / f"ibd_sensitivity_iterations_{tag}.tsv", sep="\t", index=False)
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-per-site", type=int, default=5)
    ap.add_argument("--n-iter", type=int, default=100)
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="Mantel permutations per resampling iteration")
    ap.add_argument("--n-perm-final", type=int, default=10_000,
                    help="Mantel permutations for the full-data and averaged-FST matrices")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--subregions", nargs="*", default=list(PUBLISHED))
    ap.add_argument("--florida-along-tract", action="store_true",
                    help="use along-tract rather than great-circle distance for Florida")
    args = ap.parse_args()

    if args.florida_along_tract:
        PUBLISHED["Florida"] = (PUBLISHED["Florida"][0], "along_tract_km")

    npz = np.load(NPZ_PATH, allow_pickle=True)
    G = npz["G"]
    sample_ids = pd.Series(npz["sample_ids"].astype(str)).str.strip().to_numpy()
    samples = pd.read_csv(SAMPLES_TSV, sep="\t")
    samples["affymetrix_id"] = samples["affymetrix_id"].astype(str).str.strip()

    rows = [
        run_subregion(sr, samples, G, sample_ids, args.n_per_site, args.n_iter,
                      args.n_perm, args.n_perm_final, args.seed)
        for sr in args.subregions
    ]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(rows)
    suffix = "_alongtract" if args.florida_along_tract else ""
    summary.to_csv(OUT_DIR / f"ibd_sensitivity_summary{suffix}.tsv", sep="\t", index=False)
    print("\n=== summary ===")
    print(summary.to_string(index=False))
    print(f"\nwritten to {OUT_DIR}")


if __name__ == "__main__":
    main()
