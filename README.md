# 🧬 B3GAT1 LocusZoom Benchmark

This repository contains the reproducible benchmark setup for **Figure 5B** from the paper  
> [Efficient candidate drug target discovery through proteogenomics in a Scottish cohort](https://www.nature.com/articles/s42003-025-08738-w)

using the **GCST006085 prostate cancer GWAS** dataset.

---

## Repository Contents

| File / Folder | Description |
|----------------|-------------|
| `expert_solution.R` | R-based expert solution using [`locuszoomr`](https://cran.r-project.org/web/packages/locuszoomr) |
| `solution.sh` | Wrapper script to execute the expert solution |
| `input_data/` | Contains all required local inputs (GWAS summary stats, LD data, recombination map) |
| `expert_solution.csv` | Output data table for the reproduced figure |
| `expert_solution.png` | Generated LocusZoom plot (B3GAT1 locus) |
| `problem.jsonl` | Metadata specification for benchmark ingestion |
| `Dockerfile` | Environment specification for reproducible runs |
| `paper.pdf` | Reference paper associated with this benchmark |

---

## Objective

Reproduce the **LocusZoom plot** for the *B3GAT1* locus (chr11:134.13–134.43 Mb, ±150 kb window) showing GWAS association results for prostate cancer.  
The benchmark validates the model’s ability to reconstruct locus-level association visualization from open GWAS summary statistics, LD structure, and recombination data.

---

## Analysis Overview

1. **GWAS data ingestion** from `GCST006085` (prostate cancer study).
2. **Regional filtering** for *B3GAT1* locus (±150 kb).
3. **LD integration** from pre-downloaded LDlink file (`B3GAT1_LD_r2_clean.tsv`).
4. **Recombination overlay** using UCSC `hapMapRelease24CombinedRecombMap.bw`.
5. **LocusZoom plot generation** with `locuszoomr`.

> The `r2` column in `expert_solution.csv` contains `NA` values for variants outside the *B3GAT1* ±150 kb region.  
> LD information was only computed for this local window, in line with benchmark requirements.

---

##Software Environment

| Tool | Version / Source |
|------|------------------|
| R | 4.4.0 |
| locuszoomr | CRAN release |
| data.table | ≥ 1.14 |
| ensembldb / EnsDb.Hsapiens.v75 | Bioconductor |
| rtracklayer | Bioconductor |
| ggplot2, dplyr | CRAN |
| Docker base | `rocker/r-ver:4.4.0` |

---

## Usage

### Local execution
```bash
bash solution.sh

