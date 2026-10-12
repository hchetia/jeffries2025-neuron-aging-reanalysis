# Inventory every sheet of the Supplementary Tables workbook so the analysis
# scripts can be written against the real column names.
# Run from the project root:  Rscript scripts/01_inspect_supp.R
suppressPackageStartupMessages(library(readxl))
f <- list.files("data/raw/jeffries2025", pattern = "MOESM3.*\\.xlsx$", full.names = TRUE)
if (!length(f)) {
  z <- list.files("data/raw/jeffries2025", pattern = "MOESM3.*\\.zip$", full.names = TRUE)
  stopifnot("Supplementary Tables not found in data/raw/jeffries2025" = length(z) > 0)
  unzip(z[1], exdir = "data/raw/jeffries2025/MOESM3")
  f <- list.files("data/raw/jeffries2025/MOESM3", pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)
}
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
out <- "results/tables/supp_inventory.txt"
sink(out)
for (wb in f) {
  cat("########", basename(wb), "\n")
  for (s in excel_sheets(wb)) {
    raw <- suppressMessages(read_excel(wb, sheet = s, col_names = FALSE, n_max = 8))
    full <- suppressMessages(read_excel(wb, sheet = s, col_names = FALSE))
    cat("\n=== SHEET:", s, "| rows:", nrow(full), "| cols:", ncol(full), "\n")
    print(as.data.frame(raw), right = FALSE, max = 200)
  }
}
sink()
message("Wrote ", out)
