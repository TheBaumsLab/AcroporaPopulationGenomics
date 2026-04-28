#!/usr/bin/env python3
"""
Convert a gzipped VCF to a cached NPZ genotype matrix.

Usage:
    python convert_vcf_to_npz.py --vcf path/to/input.vcf.gz --out path/to/output.npz

Requires: bcftools on PATH.
Intermediate TSV files are written to a temp directory and cleaned up automatically.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gt_io import basic_parse_report, read_gt_matrix_tsv, write_cache_npz


def _run(cmd: str) -> str:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd}\n\nSTDERR:\n{result.stderr}")
    return result.stdout.strip()


def export_vcf_to_tsv(vcf: Path, sample_names_txt: Path, gt_matrix_tsv: Path) -> dict:
    """Extract sample list and GT matrix from VCF using bcftools."""
    _run(f"bcftools query -l '{vcf}' > '{sample_names_txt}'")
    _run(f"bcftools query -f '%ID[\\t%GT]\\n' '{vcf}' > '{gt_matrix_tsv}'")
    return {
        "n_samples": int(_run(f"wc -l < '{sample_names_txt}'")),
        "n_variants": int(_run(f"wc -l < '{gt_matrix_tsv}'")),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vcf", required=True, type=Path, help="Input VCF (.vcf or .vcf.gz)")
    parser.add_argument("--out", required=True, type=Path, help="Output .npz path")
    parser.add_argument(
        "--chunksize", type=int, default=200,
        help="TSV rows per parse chunk (default: 200)",
    )
    args = parser.parse_args()

    if not args.vcf.exists():
        sys.exit(f"ERROR: VCF not found: {args.vcf}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        print(f"Extracting from {args.vcf} ...")
        stats = export_vcf_to_tsv(args.vcf, tmp / "sample_names.txt", tmp / "gt_matrix.tsv")
        print(f"  {stats['n_samples']} samples, {stats['n_variants']} variants")

        print("Parsing and writing NPZ ...")
        gt_data = read_gt_matrix_tsv(
            str(tmp / "gt_matrix.tsv"),
            str(tmp / "sample_names.txt"),
            chunksize=args.chunksize,
        )
        write_cache_npz(str(args.out), gt_data)
        print(f"Wrote: {args.out}\n")

    print(basic_parse_report(gt_data))


if __name__ == "__main__":
    main()
