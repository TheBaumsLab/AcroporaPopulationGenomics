# popgen_utils.py
"""
Shared computational utilities for A. palmata population genetics analysis.

Sections
--------
Genotype conversion
    dosage_to_genotype_array

FST
    fst_pair_with_nloci
    linearized_fst
    fst_mat_to_nm_long

Distance
    haversine_km

IBD statistics
    build_symmetric_pair_matrix
    mantel_test_pearson
    compute_ibd_stats_from_long_df

Tract ordering
    add_tract_positions
    order_sites_greedy_nn

Site-name harmonization
    assign_florida_name_family
"""
from __future__ import annotations

import re
from math import atan2, cos, radians, sin, sqrt
from typing import Optional

import allel
import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Genotype conversion
# ---------------------------------------------------------------------------

def dosage_to_genotype_array(gn: np.ndarray) -> allel.GenotypeArray:
    """
    Convert a dosage matrix (variants × samples) with values {-1, 0, 1, 2}
    into a scikit-allel GenotypeArray of shape (variants, samples, 2).

    Missing genotypes (-1) are encoded as allele value -1 in both slots.
    """
    n_variants, n_samples = gn.shape
    gt = np.full((n_variants, n_samples, 2), -1, dtype=np.int8)
    gt[gn == 0] = [0, 0]
    gt[gn == 1] = [0, 1]
    gt[gn == 2] = [1, 1]
    return allel.GenotypeArray(gt)


# ---------------------------------------------------------------------------
# FST
# ---------------------------------------------------------------------------

def fst_pair_with_nloci(
    gta: allel.GenotypeArray,
    pop1_mask: np.ndarray,
    pop2_mask: np.ndarray,
) -> tuple[float, int]:
    """
    Compute Weir-Cockerham FST between two populations defined by boolean masks.

    Returns
    -------
    fst : float
        Genome-wide FST estimate (NaN if no valid loci).
    n_loci : int
        Number of loci with finite, positive denominator.
    """
    subpops = [np.where(pop1_mask)[0], np.where(pop2_mask)[0]]
    with np.errstate(divide="ignore", invalid="ignore"):
        a, b, c = allel.weir_cockerham_fst(gta, subpops)

    num = a.sum(axis=1)
    den = (a + b + c).sum(axis=1)
    valid = np.isfinite(num) & np.isfinite(den) & (den > 0)

    if valid.sum() == 0:
        return np.nan, 0

    fst = num[valid].sum() / den[valid].sum()
    return float(fst), int(valid.sum())


def linearized_fst(fst: float) -> float:
    """Return FST / (1 - FST). Returns NaN if fst >= 1."""
    if not np.isfinite(fst) or fst >= 1:
        return np.nan
    return fst / (1.0 - fst)


def fst_mat_to_nm_long(
    fst_mat: pd.DataFrame,
    region_coords: pd.DataFrame,
    region_col: str = "region",
    lat_col: str = "lat",
    lon_col: str = "lon",
    pair_group_map: Optional[dict] = None,
) -> pd.DataFrame:
    """
    Convert a symmetric pairwise FST matrix into a long-form table with
    implied number of migrants (Nm = (1 - FST) / (4 * FST)) and geographic
    distance between region centroids.

    Parameters
    ----------
    fst_mat : pd.DataFrame
        Symmetric square matrix indexed and columned by region name.
    region_coords : pd.DataFrame
        Table with at least [region_col, lat_col, lon_col] columns.
    pair_group_map : dict, optional
        Maps region name → group label for coloring (e.g. "Northern", "MesoAm").

    Returns
    -------
    pd.DataFrame with columns: region1, region2, pair, fst, n_migrants,
        log10_n_migrants, geo_dist_km, and optionally group1/group2/pair_group.
    """
    coord_lookup = region_coords.set_index(region_col)
    rows = []

    for i, r1 in enumerate(fst_mat.index):
        for j, r2 in enumerate(fst_mat.columns):
            if i >= j:
                continue
            fst = fst_mat.loc[r1, r2]

            nm = np.nan
            log10_nm = np.nan
            if pd.notna(fst) and 0 < fst < 1:
                nm = (1 - fst) / (4 * fst)
                if nm > 0:
                    log10_nm = np.log10(nm)

            geo_dist = haversine_km(
                coord_lookup.loc[r1, lat_col], coord_lookup.loc[r1, lon_col],
                coord_lookup.loc[r2, lat_col], coord_lookup.loc[r2, lon_col],
            )

            row = {
                "region1": r1, "region2": r2,
                "pair": f"{r1} \u2013 {r2}",
                "fst": fst,
                "n_migrants": nm,
                "log10_n_migrants": log10_nm,
                "geo_dist_km": geo_dist,
            }

            if pair_group_map is not None:
                g1 = pair_group_map.get(r1, "Other")
                g2 = pair_group_map.get(r2, "Other")
                row["group1"] = g1
                row["group2"] = g2
                row["pair_group"] = g1 if g1 == g2 else "Cross-basin"

            rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Distance
# ---------------------------------------------------------------------------

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in km between two (lat, lon) points."""
    R = 6371.0
    dlat = radians(lat2 - lat1)
    dlon = radians(lon2 - lon1)
    a = (
        sin(dlat / 2) ** 2
        + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2) ** 2
    )
    return R * 2 * atan2(sqrt(a), sqrt(1 - a))


# ---------------------------------------------------------------------------
# IBD statistics
# ---------------------------------------------------------------------------

def build_symmetric_pair_matrix(
    df: pd.DataFrame,
    value_col: str,
    site1_col: str = "site1",
    site2_col: str = "site2",
) -> pd.DataFrame:
    """
    Build a symmetric square matrix from a long pairwise DataFrame.
    Diagonal is set to 0.
    """
    sites = sorted(set(df[site1_col]).union(df[site2_col]))
    mat = pd.DataFrame(np.nan, index=sites, columns=sites, dtype=float)
    for _, row in df.iterrows():
        s1, s2, val = row[site1_col], row[site2_col], row[value_col]
        mat.loc[s1, s2] = val
        mat.loc[s2, s1] = val
    np.fill_diagonal(mat.values, 0.0)
    return mat


def mantel_test_pearson(
    x_mat: pd.DataFrame,
    y_mat: pd.DataFrame,
    n_perm: int = 10_000,
    seed: int = 42,
) -> tuple[float, float]:
    """
    Mantel test using Pearson correlation on the upper triangle.

    Returns
    -------
    r_obs : float
        Observed correlation.
    p_perm : float
        Two-tailed permutation p-value.
    """
    rng = np.random.default_rng(seed)
    x = x_mat.to_numpy(dtype=float)
    y = y_mat.to_numpy(dtype=float)

    iu = np.triu_indices_from(x, k=1)
    x_vec, y_vec = x[iu], y[iu]
    valid = np.isfinite(x_vec) & np.isfinite(y_vec)
    x_vec, y_vec = x_vec[valid], y_vec[valid]

    if len(x_vec) < 2:
        return np.nan, np.nan

    r_obs = stats.pearsonr(x_vec, y_vec).statistic
    n = x.shape[0]
    perm_rs = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(n)
        y_perm_vec = y[perm][:, perm][iu][valid]
        perm_rs[i] = stats.pearsonr(x_vec, y_perm_vec).statistic

    p_perm = (np.sum(np.abs(perm_rs) >= np.abs(r_obs)) + 1) / (n_perm + 1)
    return float(r_obs), float(p_perm)


def compute_ibd_stats_from_long_df(
    df: pd.DataFrame,
    x_col: str = "geo_dist_km",
    y_col: str = "linearized_fst",
    site1_col: str = "site1",
    site2_col: str = "site2",
    n_perm: int = 10_000,
    seed: int = 42,
) -> dict:
    """
    Compute OLS (slope, intercept, R², p) and Mantel (r, p) statistics
    from a long pairwise IBD DataFrame.
    """
    x = df[x_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]

    if len(x) < 2:
        return {
            "x_col": x_col, "n_pairs": len(x),
            "slope": np.nan, "intercept": np.nan,
            "ols_r2": np.nan, "ols_p": np.nan,
            "mantel_r": np.nan, "mantel_p": np.nan,
        }

    ols = stats.linregress(x, y)
    dist_mat = build_symmetric_pair_matrix(df, x_col, site1_col=site1_col, site2_col=site2_col)
    fst_mat = build_symmetric_pair_matrix(df, y_col, site1_col=site1_col, site2_col=site2_col)
    common = dist_mat.index.intersection(fst_mat.index)
    mantel_r, mantel_p = mantel_test_pearson(
        dist_mat.loc[common, common],
        fst_mat.loc[common, common],
        n_perm=n_perm, seed=seed,
    )

    return {
        "x_col": x_col, "n_pairs": len(x),
        "slope": ols.slope, "intercept": ols.intercept,
        "ols_r2": ols.rvalue ** 2, "ols_p": ols.pvalue,
        "mantel_r": mantel_r, "mantel_p": mantel_p,
    }


# ---------------------------------------------------------------------------
# Tract ordering
# ---------------------------------------------------------------------------

def add_tract_positions(
    df: pd.DataFrame,
    site_col: str = "site",
    lat_col: str = "lat_median",
    lon_col: str = "lon_median",
) -> pd.DataFrame:
    """
    Add a cumulative along-tract distance column (``tract_km``) to an
    already-ordered site DataFrame.
    """
    out = df.copy().reset_index(drop=True)
    out["tract_km"] = 0.0
    for i in range(1, len(out)):
        prev, curr = out.iloc[i - 1], out.iloc[i]
        out.loc[i, "tract_km"] = out.loc[i - 1, "tract_km"] + haversine_km(
            prev[lat_col], prev[lon_col], curr[lat_col], curr[lon_col]
        )
    return out


def order_sites_greedy_nn(
    site_df: pd.DataFrame,
    site_col: str = "site",
    lat_col: str = "lat_median",
    lon_col: str = "lon_median",
    start_site: Optional[str] = None,
) -> pd.DataFrame:
    """
    Order sites by greedily walking from a start site to the nearest
    unvisited site (nearest-neighbor heuristic).

    Good first-pass tract ordering for curved or elbow-shaped coastlines.
    Defaults to starting at the westernmost site (minimum longitude).

    Returns the ordered DataFrame with a ``tract_km`` column appended.
    """
    df = site_df.reset_index(drop=True).copy()
    if start_site is None:
        start_idx = df[lon_col].idxmin()
    else:
        start_idx = df.index[df[site_col] == start_site][0]

    visited = [start_idx]
    remaining = set(df.index) - {start_idx}

    while remaining:
        current = df.loc[visited[-1]]
        best_idx = min(
            remaining,
            key=lambda idx: haversine_km(
                current[lat_col], current[lon_col],
                df.loc[idx, lat_col], df.loc[idx, lon_col],
            ),
        )
        visited.append(best_idx)
        remaining.remove(best_idx)

    return add_tract_positions(
        df.loc[visited].reset_index(drop=True),
        site_col=site_col, lat_col=lat_col, lon_col=lon_col,
    )


# ---------------------------------------------------------------------------
# Site-name harmonization
# ---------------------------------------------------------------------------

def assign_florida_name_family(name) -> str:
    """
    Collapse Florida reef-name variants into a single canonical site label.
    Based on regex/keyword matching of known naming conventions in the dataset.
    Returns the original name unchanged if no rule matches.
    """
    if pd.isna(name):
        return "Unknown"

    s = re.sub(r"[^a-z0-9]+", " ", str(name).strip().lower()).strip()

    rules = [
        ("sand key",       "Sand Key"),
        ("sand island",    "Sand Island"),
        ("si3",            "Sand Island"),
        ("sambo",          "Sambo"),
        ("elbow",          "Elbow"),
        ("el2",            "Elbow"),
        ("looe",           "Looe"),
        ("biscayne",       "Biscayne"),
        ("conch",          "Conch"),
        ("crf",            "Conch"),
        ("fowey",          "Fowey"),
        ("french",         "French Reef"),
        ("fr2",            "French Reef"),
        ("grecian",        "Grecian Rocks"),
        ("gr1",            "Grecian Rocks"),
        ("horseshoe",      "Horseshoe"),
        ("brew",           "Brewster"),
        ("marker3",        "Marker3"),
        ("marker 3",       "Marker3"),
        ("ball buoy",      "Ball Buoy"),
        ("ballbuoy",       "Ball Buoy"),
        ("ball paul",      "Ball Buoy"),
        ("carysfort",      "Carysfort"),
        ("cf2",            "Carysfort"),
        ("molasses",       "Molasses"),
        ("ml3",            "Molasses"),
        ("snapledge",      "Snapper Ledge"),
        ("kl4",            "Dry Rocks"),
        ("dry rocks",      "Dry Rocks"),
        ("triple a",       "Triple A"),
        ("watson",         "Watson"),
        ("trtlrks",        "Turtle Rocks"),
        ("mote",           "Mote"),
    ]
    for keyword, canonical in rules:
        if keyword in s:
            return canonical

    return str(name)
