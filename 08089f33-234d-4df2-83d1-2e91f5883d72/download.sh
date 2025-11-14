#!/usr/bin/env bash
set -euo pipefail

curl -L "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006085/29892016-GCST006085-EFO_0001663-build37.f.tsv.gz" \
  -o "./29892016-GCST006085-EFO_0001663-build37.f.tsv.gz"

