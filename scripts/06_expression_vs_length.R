# Q2/Q3. Is age-related downregulation explained by gene length or by gene function?
# Paper (Fig. 5f,g): in excitatory neurons log2(elderly/adult) rises with length in
# housekeeping (HK) genes (R2 = 0.50) but much less in neuron-specific (NS) genes (R2 = 0.20),
# read as "short housekeeping genes are selectively lost".
#
# Tests here, using the per-gene fold changes in Table S11:
#  A. reproduce the decile and gene-level fits
#  B. compare HK and NS inside the same length window (the two sets barely overlap in length)
#  C. repeat across all 15 cell types
#  D. external housekeeping lists (Eisenberg & Levanon bulk HK; Lin et al. scSEG)
#  E. drop ribosomal-protein and OXPHOS genes
# Run from project root:  Rscript scripts/06_expression_vs_length.R
suppressPackageStartupMessages(library(data.table))

fc  <- fread("data/processed/geneset_log2fc.tsv")
len <- fread("data/processed/gene_lengths.tsv")
ali <- fread("data/processed/gene_aliases.tsv")
e <- new.env(); load("data/raw/scmerge/segList.rda", envir = e)
ext_bulk <- e$segList$human$bulkRNAseqHK
ext_sc   <- e$segList$human$human_scSEG

fc[, symbol := fifelse(gene %in% len$symbol, gene, ali$symbol[match(gene, ali$alias)])]
fc <- merge(fc, len[, .(symbol, gene_length)], by = "symbol", all.x = TRUE)
fc <- fc[!is.na(gene_length) & !is.na(log2fc)]
# Drop genes whose annotated length disagrees with the paper's length decile by >5x
# (UCSC vs GENCODE gene-model differences, e.g. read-through loci); ~1.5% of genes
fc[, dec_med := as.numeric(median(gene_length)), by = .(set, length_decile, cell_type)]
message("genes dropped for length/decile disagreement: ",
        uniqueN(fc[gene_length > 5 * dec_med | gene_length < dec_med / 5, gene]))
fc <- fc[gene_length <= 5 * dec_med & gene_length >= dec_med / 5][, dec_med := NULL]
fc[, l10 := log10(gene_length)]
fc[, `:=`(in_ext_bulk = gene %in% ext_bulk | symbol %in% ext_bulk,
          in_ext_sc   = gene %in% ext_sc   | symbol %in% ext_sc)]
fc[, rp_oxphos := grepl("^(RPL|RPS)[0-9]|^(NDUF|COX[4-8]|ATP5|UQCR|SDH[ABCD])", symbol)]

fit <- function(d) {
  if (nrow(d) < 10) return(data.table(n = nrow(d), slope = NA_real_, se = NA_real_, R2 = NA_real_, p = NA_real_, rho = NA_real_))
  m <- summary(lm(log2fc ~ l10, data = d))
  data.table(n = nrow(d), slope = m$coefficients[2, 1], se = m$coefficients[2, 2], R2 = m$r.squared,
             p = m$coefficients[2, 4], rho = cor(d$log2fc, d$l10, method = "spearman"))
}
out <- list()

## A. Reproduce, excitatory neurons ("ext") -------------------------------------
ext <- fc[cell_type == "ext"]
dec_means <- ext[, .(mfc = mean(log2fc)), by = .(set, length_decile)]
A <- rbind(
  ext[, fit(.SD), by = set][, level := "gene (continuous log10 length)"],
  ext[, { m <- summary(lm(log2fc ~ length_decile)); .(n = .N, slope = coef(m)[2, 1], se = coef(m)[2, 2],
          R2 = m$r.squared, p = coef(m)[2, 4], rho = cor(log2fc, length_decile, method = "spearman")) },
      by = set][, level := "gene (decile as predictor)"],
  dec_means[, { m <- summary(lm(mfc ~ length_decile)); .(n = .N, slope = coef(m)[2, 1], se = coef(m)[2, 2],
                R2 = m$r.squared, p = coef(m)[2, 4], rho = NA_real_) }, by = set][, level := "10 decile means"])
out$A <- A[, analysis := "A_reproduce_ext"]

## B. Same length window ---------------------------------------------------------
rng <- ext[, .(lo = min(gene_length), hi = max(gene_length)), by = set]
lo <- max(rng$lo); hi <- min(rng$hi)
lo <- max(lo, ext[set == "neuron_specific", quantile(gene_length, 0.05)])
hi <- min(hi, ext[set == "housekeeping",  quantile(gene_length, 0.95)])
win <- ext[gene_length >= lo & gene_length <= hi]
B <- win[, fit(.SD), by = set][, level := sprintf("window %.0f-%.0f kb", lo / 1e3, hi / 1e3)]
anc <- summary(lm(log2fc ~ l10 * set, data = win))$coefficients
anc_add <- summary(lm(log2fc ~ l10 + set, data = win))$coefficients
message(sprintf("Length window %.0f-%.0f kb: HK n=%d, NS n=%d", lo / 1e3, hi / 1e3,
                win[set == "housekeeping", .N], win[set == "neuron_specific", .N]))
message(sprintf("  slope difference NS-HK = %.3f (p = %.3g); offset NS-HK at equal length = %.3f (p = %.3g)",
                anc["l10:setneuron_specific", 1], anc["l10:setneuron_specific", 4],
                anc_add["setneuron_specific", 1], anc_add["setneuron_specific", 4]))
out$B <- B[, analysis := "B_length_matched_ext"]
fwrite(data.table(term = rownames(anc), anc), "results/tables/expr_length_ancova_ext.tsv", sep = "\t")

# Pooled HK+NS: R2 of length alone vs set alone vs both
pooled <- ext
r2 <- data.table(model = c("length only", "set only", "length + set", "length * set"),
  R2 = c(summary(lm(log2fc ~ l10, pooled))$r.squared, summary(lm(log2fc ~ set, pooled))$r.squared,
         summary(lm(log2fc ~ l10 + set, pooled))$r.squared, summary(lm(log2fc ~ l10 * set, pooled))$r.squared))
print(r2)
fwrite(r2, "results/tables/expr_length_vs_set_R2_ext.tsv", sep = "\t")

## C. All cell types --------------------------------------------------------------
C <- fc[, fit(.SD), by = .(cell_type, set)]
Cw <- fc[gene_length >= lo & gene_length <= hi, fit(.SD), by = .(cell_type, set)]
out$C <- rbind(C[, level := "all genes"], Cw[, level := "length window"])[, analysis := "C_by_cell_type"]

## D. External housekeeping definitions (within the paper's HK set) --------------
hk <- ext[set == "housekeeping"]
message(sprintf("Paper HK genes also in Eisenberg bulk HK: %d/%d; in Lin scSEG: %d/%d",
                hk[, sum(in_ext_bulk)], nrow(hk), hk[, sum(in_ext_sc)], nrow(hk)))
D <- rbind(hk[in_ext_bulk == TRUE, fit(.SD)][, `:=`(level = "HK & Eisenberg bulk HK", set = "housekeeping")],
           hk[in_ext_bulk == FALSE, fit(.SD)][, `:=`(level = "HK, not in Eisenberg", set = "housekeeping")],
           hk[in_ext_sc == TRUE, fit(.SD)][, `:=`(level = "HK & Lin scSEG", set = "housekeeping")],
           ext[set == "neuron_specific" & !in_ext_bulk, fit(.SD)][, `:=`(level = "NS, not in Eisenberg", set = "neuron_specific")])
out$D <- D[, analysis := "D_external_HK"]

## E. Without ribosomal / OXPHOS genes -------------------------------------------
message("RP/OXPHOS genes in paper HK set: ", hk[, sum(rp_oxphos)])
E <- rbind(hk[rp_oxphos == FALSE, fit(.SD)][, level := "HK without RP/OXPHOS"],
           hk[rp_oxphos == FALSE & gene_length >= lo & gene_length <= hi, fit(.SD)][, level := "HK without RP/OXPHOS, length window"])
out$E <- E[, `:=`(set = "housekeeping", analysis = "E_drop_RP_OXPHOS")]

res <- rbindlist(out, use.names = TRUE, fill = TRUE)
setcolorder(res, c("analysis", "level", "set", "cell_type"))
fwrite(res, "results/tables/expression_vs_length_models.tsv", sep = "\t")
print(res[analysis != "C_by_cell_type"], digits = 3)
print(dcast(res[analysis == "C_by_cell_type"], cell_type ~ set + level, value.var = "slope"), digits = 2)
fwrite(fc, "data/processed/geneset_log2fc_with_length.tsv", sep = "\t")
