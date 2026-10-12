# Q2b. Do elderly and adult donors differ in technical covariates that could
# produce a length-dependent expression change (RNA quality, PMI, storage time,
# library depth/complexity, 10x chemistry batch)?
# Run from project root:  Rscript scripts/07_donor_covariates.R
suppressPackageStartupMessages(library(data.table))

d  <- fread("data/processed/donors.tsv", colClasses = list(character = "case_id"))
qc <- fread("data/processed/library_qc.tsv")
setnames(qc, c("library", "n_cells", "mean_reads", "mean_genes", "mean_mt_pct", "batch"))
qc[, case_id := sub("_.*", "", library)]
qc[case_id == "0950", case_id := "09-50"]
# S3 QC already lists only the libraries kept after the authors' batch exclusions
lib <- qc[, .(n_libraries = .N, cells = sum(n_cells),
              mean_reads = weighted.mean(mean_reads, n_cells),
              mean_genes = weighted.mean(mean_genes, n_cells),
              mean_mt_pct = weighted.mean(mean_mt_pct, n_cells)), by = case_id]
d <- merge(d, lib, by = "case_id", all.x = TRUE)
d[, genes_per_read := mean_genes / mean_reads]
fwrite(d, "results/tables/donor_covariates.tsv", sep = "\t")

ae <- d[group %in% c("adult", "elderly")]
vars <- c("rin", "pmi_h", "storage_yrs", "mean_reads", "mean_genes", "genes_per_read", "mean_mt_pct", "n_libraries")
tab <- rbindlist(lapply(vars, function(v) {
  x <- ae[group == "adult", get(v)]; y <- ae[group == "elderly", get(v)]
  data.table(variable = v, adult_median = median(x, na.rm = TRUE), elderly_median = median(y, na.rm = TRUE),
             n_adult = sum(!is.na(x)), n_elderly = sum(!is.na(y)),
             wilcox_p = suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value),
             spearman_rho_with_age = suppressWarnings(cor(ae$age, ae[[v]], method = "spearman", use = "complete.obs")),
             spearman_p = suppressWarnings(cor.test(ae$age, ae[[v]], method = "spearman", exact = FALSE)$p.value))
}))
print(tab, digits = 3)
fwrite(tab, "results/tables/donor_covariates_by_group.tsv", sep = "\t")

# Batch composition: libraries per sequencing batch, by group
bat <- merge(qc, d[, .(case_id, group)], by = "case_id")[group %in% c("adult", "elderly")]
bt <- dcast(bat[, .N, by = .(batch, group)], batch ~ group, value.var = "N", fill = 0)
print(bt)
fwrite(bt, "results/tables/library_batches_by_group.tsv", sep = "\t")
