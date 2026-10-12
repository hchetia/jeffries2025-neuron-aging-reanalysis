#!/usr/bin/env bash
# Download public supplementary data for Jeffries et al. 2025 (Nature 646:657).
# Run from the project root:  bash scripts/00_download.sh
set -euo pipefail
OUT="data/raw/jeffries2025"; mkdir -p "$OUT"
BASE="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09435-8/MediaObjects/41586_2025_9435"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)"

get() {  # get <MOESM id> <candidate extensions...>
  local id=$1; shift
  for ext in "$@"; do
    local f="$OUT/41586_2025_9435_${id}_ESM.${ext}"
    if curl -fsSL -A "$UA" -o "$f" "${BASE}_${id}_ESM.${ext}"; then
      echo "OK   $f ($(du -h "$f" | cut -f1))"; return 0
    fi
    rm -f "$f"
  done
  echo "FAIL $id - download manually from the article's Supplementary information section" >&2
}

get MOESM1 pdf          # Supplementary Information (Supp. Figs, methods)
get MOESM3 xlsx zip     # Supplementary Tables 1-16
# Source data / extended data, if present
for i in 4 5 6 7 8 9 10 11 12 13 14 15; do get "MOESM$i" xlsx zip csv >/dev/null 2>&1 || true; done
ls -lh "$OUT"

# Annotation used by the reanalysis (both small, from GitHub)
mkdir -p data/raw/annotables data/raw/scmerge
# Ensembl GRCh37 gene coordinates (S13 SNV calls are on hg19/GRCh37)
[ -s data/raw/annotables/grch37.rda ] || curl -fsSL -o data/raw/annotables/grch37.rda \
  https://raw.githubusercontent.com/stephenturner/annotables/master/data/grch37.rda
# External housekeeping lists: Eisenberg & Levanon 2013 bulk HK; Lin et al. 2019 scSEG
[ -s data/raw/scmerge/segList.rda ] || curl -fsSL -o data/raw/scmerge/segList.rda \
  https://raw.githubusercontent.com/SydneyBioX/scMerge/master/data/segList.rda
ls -lh data/raw/annotables data/raw/scmerge
