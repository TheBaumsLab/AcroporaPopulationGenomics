# AcroporaPopulationGenomics

Code and data accompanying:

> **Population structure and gene flow in the endangered Caribbean reef-building coral, *Acropora palmata***  
> Iliana B. Baums, Nicolas S. Locatelli, Kim L. de Luca, Sheila A. Kitchen  

Raw genotype calls and full sample metadata are available through the STAGdb Galaxy resource at [coralsnp.uol.de](https://coralsnp.uol.de). Datasets are available in gzip-compressed VCF format at [coralsnp.uol.de/galaxy/libraries](https://coralsnp.uol.de/galaxy/libraries) under "Acropora palmata popgen datasets".

---

## Repository structure

```
├── filtered_calls/          # Filtered VCF files (datasets 4 & 6, see below)
├── popgen_notebooks/        # Analysis notebooks
├── scripts/                 # Shared Python utilities
├── slurm/                   # SLURM submission scripts (ADMIXTURE, STRUCTURE, liftover)
├── resources/               # Auxiliary data (liftover chain file, IHO shapefiles)
└── requirements.txt         # Python dependencies
```

---

## Datasets

The manuscript describes six datasets produced by progressive filtering of the raw STAGdb VCF. The two VCF files in `filtered_calls/` correspond to:

| File | Dataset | Description |
|---|---|---|
| `allsamples_lifted_to_jaAcrPala1.3_conffiltered_indvmissfiltered_downsampled_maffiltered_ldfiltered.vcf.gz` | Dataset 4 | Fully filtered: probes, CONF, site/sample missingness, LD-pruned, MAF > 0.05, geographically downsampled (≤ 60 genets/region). |
| `allsamples_lifted_to_jaAcrPala1.3_noconf_indvmissfiltered_maffiltered.vcf.gz` | Dataset 6 | No CONF filter. Otherwise same as dataset 4. |

---

## Notebooks

All notebooks are in `popgen_notebooks/`.

### SNP chip analysis

| Notebook | Analysis |
|---|---|
| `1_microarray_liftover_to_jaAcrPala1.3.ipynb` | Genome coordinate liftover from *A. digitifera* to *A. palmata* (GCF_964030605.1) using Progressive Cactus and GATK LiftoverVcf |
| `2_snpchip_filtering_and_PCA.ipynb` | SNP and sample filtering pipeline; produces datasets 1–6; principal component analysis (Figure 1, Figures S1–S2) |
| `3_admixture_and_structure_snpchip_data.ipynb` | STRUCTURE and ADMIXTURE analyses; determination of optimal K; bar plot on map (Figure 2A, Figures S3–S4) |
| `4_feems_snpchip_data.ipynb` | Effective migration surface analysis with FEEMS (Figure 2B) |
| `5_fst_table.ipynb` | Pairwise FST visualization |
| `5a_fst_calculations.R` | Bootstrapped pairwise FST between regions using StAMPP (Table 2) |
| `6_KING_plots.ipynb` | Within-site kinship distributions from PLINK2 KING estimates (Figure 4) |
| `7_inbreeding_and_het.ipynb` | Regional observed heterozygosity and inbreeding coefficients (FIS) from VCFtools (Table 3) |
|---|---|
| `01_sample_summary.ipynb` | Sample inventory by region; flags facility/nursery samples (`is_facility`); region-level summary table |
| `02_fst_migrants.ipynb` | Pairwise FST between regions; migration rate estimates under Wright's island model (Nm) Table S2 |
| `03_ibd_analysis.ipynb` | Isolation-by-distance along the Greater Antilles and Florida Reef Tract; site merging, regional MAF filtering, OLS regression, Mantel test, combined IBD figure (Figure 3); writes `apal_ibd_samples.tsv`, `ibd_fst_greater_antilles.tsv`, `ibd_fst_florida.tsv` |

Notebooks 01–03 require setting the following environment variables before launching Jupyter:

```bash
export APAL_PROJECT_DIR=/path/to/your/project   # directory containing table_S1_metadata.tsv
export APAL_DATA_DIR=/path/to/your/data         # directory containing the .npz genotype cache
```

---

## Scripts

| Script | Description |
|---|---|
| `scripts/convert_vcf_to_npz.py` | Converts a filtered VCF to a compressed int8 dosage matrix (`.npz`) used by the IBD notebooks. Requires `bcftools ≥ 1.17`. |
| `scripts/gt_io.py` | Genotype I/O utilities: loading/saving the `.npz` cache, sample/variant subsetting. |
| `scripts/popgen_utils.py` | Shared computational utilities: Weir–Cockerham FST, linearized FST, haversine distance, Mantel test (`skbio`), IBD statistics, tract ordering, site-name harmonization. |

---

## Requirements

Python dependencies are listed in `requirements.txt`. Install with:

```bash
pip install -r requirements.txt
```

One additional tool is required but not pip-installable:

- **bcftools ≥ 1.17** — used by `convert_vcf_to_npz.py`
  - conda: `conda install -c bioconda bcftools`
  - brew: `brew install bcftools`

---

## Citations

> Baums IB _et al_. (2026) Population structure and gene flow in the endangered Caribbean reef-building coral, *Acropora palmata*. https://doi.org/10.64898/2026.04.15.718759  
> Kitchen SA _et al_. (2020) STAGdb: a 30K SNP genotyping array and interactive database for *Acropora* corals and their dinoflagellate symbionts. *Scientific Reports* 10, 12488. https://doi.org/10.1038/s41598-020-69101-z
