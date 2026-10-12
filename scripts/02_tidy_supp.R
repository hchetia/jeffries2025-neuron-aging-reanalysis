# Convert the Jeffries et al. 2025 supplementary workbooks into tidy TSVs.
# Run from the project root:  Rscript scripts/02_tidy_supp.R
suppressPackageStartupMessages({ library(readxl); library(data.table) })

find_xlsx <- function(name) {
  f <- list.files("data/raw/jeffries2025", pattern = paste0("^", name), full.names = TRUE, recursive = TRUE)
  if (!length(f)) stop("Cannot find ", name, " under data/raw/jeffries2025")
  f[1]
}
out <- function(x, f) { fwrite(x, file.path("data/processed", f), sep = "\t"); message("wrote ", f, " (", nrow(x), " rows)") }
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
num <- function(x) suppressWarnings(as.numeric(gsub("[,†‡]", "", x)))

## S1: donors -----------------------------------------------------------------
s1 <- as.data.table(read_excel(find_xlsx("S1_"), sheet = 1, col_types = "text"))
setnames(s1, c("case_id", "sex", "age", "rna_nuclei", "dna_nuclei", "merfish_cells",
               "pmi_h", "cause_of_death", "rin", "storage_yrs", "brain_bank"))
s1 <- s1[!is.na(case_id) & !is.na(num(age))]
for (v in c("age", "rna_nuclei", "dna_nuclei", "pmi_h", "rin", "storage_yrs")) set(s1, j = v, value = num(s1[[v]]))
# Group definitions used in the paper's two-group DE (S6 legend): adult 10<age<=60, elderly >60
s1[, group := fifelse(age < 10, "infant", fifelse(age <= 60, "adult", "elderly"))]
out(s1, "donors.tsv")

## S11: housekeeping / neuron-specific gene fold changes ----------------------
read_geneset <- function(sheet, set) {
  raw <- read_excel(find_xlsx("S11_"), sheet = sheet, col_names = FALSE, col_types = "text")
  ct <- unlist(raw[2, -(1:2)])
  d <- as.data.table(raw[-(1:2), ])
  setnames(d, c("gene", "length_decile", ct))
  d <- melt(d, id.vars = c("gene", "length_decile"), variable.name = "cell_type", value.name = "log2fc")
  d[, `:=`(length_decile = as.integer(length_decile), log2fc = num(log2fc), set = set)]
  d[!is.na(gene)]
}
gs <- rbind(read_geneset("Housekeeping Genes", "housekeeping"),
            read_geneset("Neuron-Specific Genes", "neuron_specific"))
out(gs, "geneset_log2fc.tsv")

## S13: per-neuron summary and per-SNV calls ----------------------------------
summ <- as.data.table(read_excel(find_xlsx("S13_"), sheet = "Summary", skip = 3, col_types = "text"))
setnames(summ, c("sample", "donor", "age", "n_ssnv_called", "est_ssnv_genome", "sensitivity",
                 "est_A1", "est_A2", "frac_A1", "frac_A2"))
summ <- summ[!is.na(sample) & !is.na(num(age))]
for (v in setdiff(names(summ), c("sample", "donor"))) set(summ, j = v, value = num(summ[[v]]))
out(summ, "neuron_summary.tsv")

snv <- as.data.table(read_excel(find_xlsx("S13_"), sheet = "Detailed information", skip = 2, col_types = "text"))
setnames(snv, make.names(names(snv)))
snv[, `:=`(Start = as.integer(Start), End = as.integer(End))]
out(snv, "ssnv_calls.tsv")

## S16: published decile-level SNV rates --------------------------------------
raw <- read_excel(find_xlsx("S16_"), sheet = 1, col_names = FALSE, col_types = "text")
hdr <- which(raw[[1]] == "Decile")
lab <- raw[[1]][hdr - 1]
s16 <- rbindlist(lapply(seq_along(hdr), function(i) {
  rows <- (hdr[i] + 1):(hdr[i] + 10)
  data.table(set = lab[i], decile = as.integer(raw[[1]][rows]), total_bp = num(raw[[2]][rows]),
             mutations = num(raw[[3]][rows]), neurons = num(raw[[4]][rows]), rate = num(raw[[5]][rows]))
}))
out(s16, "s16_decile_rates.tsv")

## S6: elderly vs adult DEGs (significant only) -------------------------------
s6 <- as.data.table(read_excel(find_xlsx("S6_"), sheet = "DEGs", skip = 3, col_types = "text"))
setnames(s6, c("gene", "log2fc", "pct_elderly", "pct_adult", "p", "q", "cell_type"))
s6 <- s6[!is.na(gene)]
for (v in c("log2fc", "pct_elderly", "pct_adult", "p", "q")) set(s6, j = v, value = num(s6[[v]]))
out(s6, "degs_elderly_vs_adult.tsv")

## S3: per-library QC ---------------------------------------------------------
qc <- as.data.table(read_excel(find_xlsx("S3_"), sheet = "QC", col_types = "text"))
out(qc, "library_qc.tsv")
