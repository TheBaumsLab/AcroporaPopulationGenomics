# gt_io.py
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import pandas as pd


GT_TO_INT8 = {
    "0/0": 0,
    "0/1": 1,
    "1/0": 1,
    "1/1": 2,
    "./.": -1,
    "0|0": 0,
    "0|1": 1,
    "1|0": 1,
    "1|1": 2,
    ".|.": -1,
}


@dataclass
class GTData:
    ax_ids: np.ndarray        # shape (n_snps,), dtype=object/str
    sample_ids: np.ndarray    # shape (n_samples,), dtype=object/str
    G: np.ndarray             # shape (n_snps, n_samples), dtype=int8, values {0,1,2,-1}


def read_sample_ids(sample_names_path: str) -> np.ndarray:
    """Read one sample ID per line; strips whitespace."""
    if not os.path.exists(sample_names_path):
        raise FileNotFoundError(f"Sample names file not found: {sample_names_path}")
    with open(sample_names_path, "r") as f:
        sample_ids = [line.strip() for line in f if line.strip() != ""]
    if len(sample_ids) == 0:
        raise ValueError(f"No sample IDs read from {sample_names_path}")
    return np.array(sample_ids, dtype=object)


def _validate_tokens(df: pd.DataFrame) -> List[str]:
    """Return sorted list of unique tokens found in genotype cells (excluding AX ID col)."""
    tokens = pd.unique(df.iloc[:, 1:].values.ravel())
    tokens = [t for t in tokens if isinstance(t, str)]
    return sorted(set(tokens))


def read_gt_matrix_tsv(
    gt_path: str,
    sample_names_path: str,
    *,
    expect_n_samples: Optional[int] = None,
    chunksize: Optional[int] = 200,
) -> GTData:
    """
    Read a no-header TSV:
      col0 = AX ID
      col1.. = GT tokens (0/0, 0/1, 1/1, ./.)

    Uses chunked parsing to stay memory-friendly.
    """
    if not os.path.exists(gt_path):
        raise FileNotFoundError(f"GT matrix not found: {gt_path}")

    sample_ids = read_sample_ids(sample_names_path)

    # Peek first 200 rows to validate column count and token set
    peek = pd.read_csv(gt_path, sep="\t", header=None, dtype=str, nrows=200)
    n_samples_in_matrix = peek.shape[1] - 1

    if expect_n_samples is not None and n_samples_in_matrix != expect_n_samples:
        raise ValueError(f"Matrix has {n_samples_in_matrix} samples, expected {expect_n_samples}")

    if len(sample_ids) != n_samples_in_matrix:
        raise ValueError(
            f"Sample IDs length ({len(sample_ids)}) does not match matrix columns ({n_samples_in_matrix}).\n"
            f"Check that {sample_names_path} is complete and ordered correctly."
        )

    tokens = _validate_tokens(peek)
    bad = [t for t in tokens if t not in GT_TO_INT8]
    if bad:
        raise ValueError(f"Unexpected genotype tokens found (first 200 rows): {bad}\nAll tokens: {tokens}")

    ax_list: List[str] = []
    G_list: List[np.ndarray] = []

    reader = pd.read_csv(gt_path, sep="\t", header=None, dtype=str, chunksize=chunksize)
    for chunk in reader:
        ax_list.extend(chunk.iloc[:, 0].astype(str).tolist())

        gt = chunk.iloc[:, 1:].to_numpy(dtype=str, copy=False)
        gt_num = np.empty(gt.shape, dtype=np.int8)
        gt_num[:] = -128  # sentinel for unmapped

        for tok, code in GT_TO_INT8.items():
            gt_num[gt == tok] = np.int8(code)

        if (gt_num == -128).any():
            bad = np.unique(gt[gt_num == -128]).tolist()
            raise ValueError(f"Unexpected genotype tokens in chunk: {bad}")

        G_list.append(gt_num)

    ax_ids = np.array(ax_list, dtype=object)
    G = np.vstack(G_list).astype(np.int8, copy=False)

    return GTData(ax_ids=ax_ids, sample_ids=sample_ids, G=G)


def write_cache_npz(out_npz: str, data: GTData) -> None:
    """Save compact cache: AX IDs + sample IDs + genotype matrix."""
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    np.savez_compressed(
        out_npz,
        ax_ids=data.ax_ids,
        sample_ids=data.sample_ids,
        G=data.G,
    )


def load_cache_npz(npz_path: str) -> GTData:
    """Load a cached NPZ. Accepts either 'ax_ids' or 'variant_ids' as the SNP key."""
    if not os.path.exists(npz_path):
        raise FileNotFoundError(npz_path)
    z = np.load(npz_path, allow_pickle=True)
    if "ax_ids" in z:
        ax_ids = z["ax_ids"]
    elif "variant_ids" in z:
        ax_ids = z["variant_ids"]
    else:
        raise KeyError(f"No variant ID field found in {npz_path}")
    return GTData(
        ax_ids=ax_ids,
        sample_ids=z["sample_ids"],
        G=z["G"].astype(np.int8, copy=False),
    )


def basic_parse_report(data: GTData) -> str:
    """Return a human-readable summary of a GTData object."""
    n_snps, n_samples = data.G.shape
    missing_rate = float((data.G == -1).mean())

    ax_series = pd.Series(data.ax_ids)
    dup_mask = ax_series.duplicated(keep=False)
    n_dups = int(dup_mask.sum())
    dup_unique = int(ax_series[dup_mask].nunique()) if n_dups else 0

    vals, counts = np.unique(data.G, return_counts=True)
    counts_map = {int(v): int(c) for v, c in zip(vals, counts)}
    tok_map = {-1: "./.", 0: "0/0", 1: "0/1", 2: "1/1"}
    tok_counts = {tok_map.get(k, str(k)): v for k, v in counts_map.items()}

    return "\n".join([
        f"n_snps:                {n_snps}",
        f"n_samples:             {n_samples}",
        f"overall_missing_rate:  {missing_rate:.6f}",
        f"AX_duplicates_cells:   {n_dups} (unique duplicated AX IDs: {dup_unique})",
        f"genotype_counts:       {tok_counts}",
    ])
