#!/usr/bin/env python3
"""
Pooled ("all sites together") isolation-by-distance analysis.

Reviewer 2 (L. 171-173) asks why IBD was evaluated separately along the Greater
Antilles and along the Florida Reef Tract rather than over all sites at once.
This script runs the pooled analysis so that the answer can be supported with
numbers rather than asserted.

It pools the same curated sites used in the manuscript (Greater Antilles + Florida,
wild genets only, >5 genets per site, 10-km site merging), computes pairwise WC84
FST on a common locus panel, and then:

  1. Mantel test / OLS on all site pairs pooled.
  2. The same statistics computed separately for within-Greater-Antilles pairs,
     within-Florida pairs, and between-subregion pairs.
  3. A partial Mantel test of geographic distance vs. linearized FST controlling for
     a binary "different subregion" design matrix -- i.e. does distance still explain
     genetic differentiation once the discrete east-west split is accounted for?

The expected and reported pattern is that the pooled Mantel correlation is dominated
by the between-subregion contrast (a step, not a slope), which is why the manuscript
tests IBD within each linear habitat instead.

Inputs (paths follow 03_ibd_analysis.ipynb; set the env vars below):
  APAL_PROJECT_DIR  folder with processed_data_03_apal_ibd_samples.tsv and
                    processed_data_03a_apal_ibd_fst_greater_antilles.tsv
  APAL_DATA_DIR     folder with mac_ld.gt_matrix.int8.npz
                    (generate it from dataset 3 with scripts/convert_vcf_to_npz.py)
Output: $APAL_OUT_DIR/pooled  (default <project>/output/ibd_robustness/pooled)
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
OUT_DIR = RESULTS_DIR / "pooled"

MAF_THRESH = 0.05  # regional MAF, matching REGIONAL_MAF_MIN in the published notebook
EARTH_R_KM = 6371.0


# --- shared FST / Mantel helpers (kept in-file so the script is self-contained) ---


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
    """Per-subregion locus filter (regional MAF >= 0.05 over called alleles), exactly
    as in cells 34/61 of 03_ibd_analysis.ipynb."""
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


# --- pooled-analysis specifics ---


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_R_KM * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def partial_mantel(x_mat, y_mat, z_mat, n_perm=10_000, rng=None):
    """
    Partial Mantel (Smouse et al. 1986): correlation of x and y after regressing both
    on the control matrix z, with significance from permuting the rows/columns of y.
    """
    rng = np.random.default_rng(42) if rng is None else rng
    x, y, z = (np.asarray(m, dtype=float) for m in (x_mat, y_mat, z_mat))
    iu = np.triu_indices_from(x, k=1)

    def resid(v, w):
        keep = np.isfinite(v) & np.isfinite(w)
        b = np.polyfit(w[keep], v[keep], 1)
        out = np.full_like(v, np.nan)
        out[keep] = v[keep] - (b[0] * w[keep] + b[1])
        return out

    xv, yv, zv = x[iu], y[iu], z[iu]
    valid = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(zv)
    rx = resid(xv, zv)[valid]
    r_obs = stats.pearsonr(rx, resid(yv, zv)[valid]).statistic

    n = x.shape[0]
    perm_rs = np.empty(n_perm)
    for i in range(n_perm):
        p = rng.permutation(n)
        yp = y[p][:, p][iu]
        perm_rs[i] = stats.pearsonr(rx, resid(yp, zv)[valid]).statistic
    p_val = (np.sum(np.abs(perm_rs) >= abs(r_obs)) + 1) / (n_perm + 1)
    return r_obs, p_val


def ols_block(df, label):
    x = df["geo_dist_km"].to_numpy(float)
    y = df["linearized_fst"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    r = stats.linregress(x[ok], y[ok])
    return {
        "stratum": label,
        "n_pairs": int(ok.sum()),
        "geo_dist_km_min": float(np.nanmin(x[ok])),
        "geo_dist_km_max": float(np.nanmax(x[ok])),
        "fst_mean": float(np.nanmean(df["fst"].to_numpy(float)[ok])),
        "ols_slope": r.slope,
        "ols_r2": r.rvalue ** 2,
        "ols_p": r.pvalue,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    npz = np.load(NPZ_PATH, allow_pickle=True)
    G = npz["G"]
    sample_ids = pd.Series(npz["sample_ids"].astype(str)).str.strip().to_numpy()

    s = pd.read_csv(SAMPLES_TSV, sep="\t")
    s["affymetrix_id"] = s["affymetrix_id"].astype(str).str.strip()
    s = s.drop_duplicates(["mlg", "site_merged"]).reset_index(drop=True)

    id_to_col = pd.Series(np.arange(len(sample_ids)), index=sample_ids)
    s["col"] = s["affymetrix_id"].map(id_to_col).astype(int)

    cols = s["col"].to_numpy()
    mask = locus_mask_for(G[:, cols])
    gta = dosage_to_genotype_array(G[mask][:, cols].astype(np.int16))
    print(f"pooled panel: {len(s)} genets, {int(mask.sum())}/{G.shape[0]} loci")

    site_coords = (
        s.groupby("site_merged")
        .agg(subregion=("subregion_ibd", "first"), region=("region", "first"),
             lat=("lat_clean", "median"), lon=("lon_clean", "median"),
             n_genets=("mlg", "nunique"))
        .reset_index()
        .sort_values("site_merged")
        .reset_index(drop=True)
    )
    sites = site_coords["site_merged"].tolist()
    site_rows = {x: np.where(s["site_merged"].to_numpy() == x)[0] for x in sites}
    print(f"pooled panel: {len(sites)} sites "
          f"({(site_coords['subregion'] == 'Florida').sum()} Florida, "
          f"{(site_coords['subregion'] == 'GreaterAntilles').sum()} Greater Antilles)")

    rows = []
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            s1, s2 = sites[i], sites[j]
            fst, nl = fst_pair_with_nloci(gta, site_rows[s1], site_rows[s2])
            a, b = site_coords.iloc[i], site_coords.iloc[j]
            rows.append({
                "site1": s1, "site2": s2,
                "subregion1": a["subregion"], "subregion2": b["subregion"],
                "region1": a["region"], "region2": b["region"],
                "n_genets_site1": int(a["n_genets"]), "n_genets_site2": int(b["n_genets"]),
                "fst": fst, "n_loci": nl,
                "geo_dist_km": haversine_km(a["lat"], a["lon"], b["lat"], b["lon"]),
            })
    pairs = pd.DataFrame(rows)
    pairs["linearized_fst"] = pairs["fst"] / (1 - pairs["fst"])
    pairs["same_subregion"] = pairs["subregion1"] == pairs["subregion2"]
    pairs["stratum"] = np.where(
        ~pairs["same_subregion"], "between_subregion", "within_" + pairs["subregion1"]
    )

    # ---- validation against the published Greater Antilles distances
    pub = pd.read_csv(FST_GA_TSV, sep="\t")
    chk = pub.merge(pairs, on=["site1", "site2"], suffixes=("_pub", ""))
    if len(chk):
        print(f"geo-distance check vs published GA table: "
              f"max |delta| = {(chk['geo_dist_km_pub'] - chk['geo_dist_km']).abs().max():.3f} km")

    # ---- pooled and stratified statistics
    geo = square_from_long(pairs, "geo_dist_km", sites)
    fstm = square_from_long(pairs, "linearized_fst", sites)
    sub_design = square_from_long(
        pairs.assign(diff_sub=(~pairs["same_subregion"]).astype(float)), "diff_sub", sites
    )

    m_r, m_p = mantel_test_pearson(geo, fstm, n_perm=args.n_perm, rng=rng)
    pm_r, pm_p = partial_mantel(geo, fstm, sub_design, n_perm=args.n_perm, rng=rng)
    d_r, d_p = mantel_test_pearson(sub_design, fstm, n_perm=args.n_perm, rng=rng)

    blocks = [ols_block(pairs, "pooled_all_sites")]
    for lab, sub in pairs.groupby("stratum"):
        blocks.append(ols_block(sub, lab))
    blocks_df = pd.DataFrame(blocks)

    summary = pd.DataFrame([{
        "n_sites": len(sites),
        "n_pairs": len(pairs),
        "n_genets": len(s),
        "n_loci": int(mask.sum()),
        "mantel_r_pooled": m_r,
        "mantel_p_pooled": m_p,
        "mantel_r_subregion_split_only": d_r,
        "mantel_p_subregion_split_only": d_p,
        "partial_mantel_r_geo_given_split": pm_r,
        "partial_mantel_p_geo_given_split": pm_p,
    }])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT_DIR / "ibd_pooled_pairs.tsv", sep="\t", index=False)
    blocks_df.to_csv(OUT_DIR / "ibd_pooled_strata.tsv", sep="\t", index=False)
    summary.to_csv(OUT_DIR / "ibd_pooled_summary.tsv", sep="\t", index=False)

    print("\n=== pooled Mantel ===")
    print(summary.to_string(index=False))
    print("\n=== stratified OLS ===")
    print(blocks_df.to_string(index=False))
    print(f"\nwritten to {OUT_DIR}")


if __name__ == "__main__":
    main()
