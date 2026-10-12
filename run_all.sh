#!/usr/bin/env bash
# Full pipeline. Run from the project root after `Rscript scripts/install_packages.R`.
set -euo pipefail
[ -d data/raw/jeffries2025 ] && ls data/raw/jeffries2025 | grep -q . || bash scripts/00_download.sh
if [ -f data/raw/jeffries2025/41586_2025_9435_MOESM3_ESM.zip ] && [ ! -d data/raw/jeffries2025/MOESM3 ]; then
  unzip -q data/raw/jeffries2025/41586_2025_9435_MOESM3_ESM.zip -d data/raw/jeffries2025/MOESM3
fi
mkdir -p data/processed results/tables results/figures
for s in 02_tidy_supp 03_gene_lengths 04_reproduce_s16 05_mutation_vs_length 06_expression_vs_length 07_donor_covariates 08_figures; do
  echo "== $s"; Rscript "scripts/$s.R"
done
