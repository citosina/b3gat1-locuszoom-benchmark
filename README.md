# B3GAT1 LocusZoom Benchmark

This repo contains the reproducible benchmark setup for reproducing **Figure 5B** from the paper https://www.nature.com/articles/s42003-025-08738-w with the data source being the GCST006085 prostate cancer GWAS study.

## Contents
- `expert_solution.R`: R-based expert solution using `locuszoomr`
- `solution.sh`: wrapper script for execution
- `input_data/`: contains all required local inputs (GWAS, LD, recombination)
- `expert_solution.csv` and `expert_solution.png`: generated outputs
- `problem.jsonl`: metadata for benchmark ingestion
- `Dockerfile`: environment specification for reproducible runs

## Usage
Run the solution locally:
```bash
bash solution.sh



Ensure the following files exist in input_data/:

29892016-GCST006085-EFO_0001663-build37.f.tsv.gz

B3GAT1_LD_r2_clean.tsv

recomb_hg19_hapMapRelease24CombinedRecombMap.bw

The GWAS summary file is too large for GitHub and must be manually placed in input_data/.

Outputs

expert_solution.csv: tabular output for benchmark ingestion

expert_solution.png: LocusZoom plot reproducing the target figure

